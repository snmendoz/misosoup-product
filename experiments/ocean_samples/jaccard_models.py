#!/usr/bin/env python3
"""Pairwise Jaccard distances among gapseq metabolic models.

For each model set, Jaccard distance is computed on the set of SBML reaction IDs:

    d_J(A, B) = 1 - |A ∩ B| / |A ∪ B|

Only the strict upper triangle is evaluated. The symmetric lower triangle is
filled from the computed values, and the diagonal is zero.

The script is designed for the TARA Chile 1,375-MAG model collections and
validates that the same model IDs occur in all three sets before computing any
matrix.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np


DEFAULT_SETS = (
    ("draft", "Models/gapseq/models/draft"),
    ("complete_medium", "Models/gapseq/models/gapfilled/complete_medium"),
    ("LS2N_medium", "Models/gapseq/models/gapfilled/LS2N_medium"),
)


def natural_key(value: str):
    return tuple(
        int(token) if token.isdigit() else token.lower()
        for token in re.split(r"(\d+)", value)
    )


def model_key(path: Path) -> str:
    """Return a stable MAG/model key from a model filename."""
    stem = path.stem.strip()
    suffixes = (
        "-draft", "_draft", ".draft",
        "-gapfilled", "_gapfilled", ".gapfilled",
    )
    changed = True
    while changed:
        changed = False
        lower = stem.lower()
        for suffix in suffixes:
            if lower.endswith(suffix):
                stem = stem[: -len(suffix)]
                changed = True
                break
    return stem


def discover_models(directory: Path, expected_models: int) -> dict[str, Path]:
    directory = directory.expanduser().resolve()
    if not directory.is_dir():
        raise FileNotFoundError(f"Model directory not found: {directory}")

    paths = sorted(
        [
            p
            for p in directory.iterdir()
            if p.is_file() and p.suffix.lower() in {".xml", ".sbml"}
        ],
        key=lambda p: natural_key(p.name),
    )

    if len(paths) != expected_models:
        raise ValueError(
            f"{directory}: expected {expected_models} SBML models, found {len(paths)}"
        )

    mapping: dict[str, Path] = {}
    duplicates: dict[str, list[Path]] = {}
    for path in paths:
        key = model_key(path)
        if key in mapping:
            duplicates.setdefault(key, [mapping[key]]).append(path)
        else:
            mapping[key] = path

    if duplicates:
        detail = "; ".join(
            f"{key}: {[p.name for p in files]}"
            for key, files in list(duplicates.items())[:10]
        )
        raise ValueError(f"Duplicate normalized model IDs in {directory}: {detail}")

    return mapping


def validate_model_sets(
    data_root: Path,
    expected_models: int,
    set_specs=DEFAULT_SETS,
) -> tuple[list[str], dict[str, dict[str, Path]]]:
    mappings: dict[str, dict[str, Path]] = {}
    for set_name, relative_path in set_specs:
        mappings[set_name] = discover_models(
            data_root / relative_path,
            expected_models=expected_models,
        )

    reference_name = set_specs[0][0]
    reference_ids = set(mappings[reference_name])

    for set_name, _ in set_specs[1:]:
        ids = set(mappings[set_name])
        missing = sorted(reference_ids - ids, key=natural_key)
        extra = sorted(ids - reference_ids, key=natural_key)
        if missing or extra:
            raise ValueError(
                f"Model IDs do not match between {reference_name} and {set_name}. "
                f"Missing in {set_name}: {missing[:20]}; "
                f"extra in {set_name}: {extra[:20]}"
            )

    ordered_ids = sorted(reference_ids, key=natural_key)
    return ordered_ids, mappings


def read_reaction_ids(path: Path) -> set[str]:
    """Stream an SBML file and return its reaction IDs."""
    reaction_ids: set[str] = set()
    try:
        for _event, elem in ET.iterparse(path, events=("end",)):
            if elem.tag.rsplit("}", 1)[-1] == "reaction":
                reaction_id = elem.attrib.get("id")
                if reaction_id:
                    reaction_ids.add(reaction_id)
            elem.clear()
    except ET.ParseError as exc:
        raise ValueError(f"Invalid XML/SBML file {path}: {exc}") from exc

    if not reaction_ids:
        raise ValueError(f"No SBML reactions found in {path}")
    return reaction_ids


def reaction_sets_to_bitmasks(
    reaction_sets: list[set[str]],
) -> tuple[list[int], int]:
    universe = sorted(set().union(*reaction_sets))
    bit_index = {reaction_id: idx for idx, reaction_id in enumerate(universe)}

    masks: list[int] = []
    for reactions in reaction_sets:
        mask = 0
        for reaction_id in reactions:
            mask |= 1 << bit_index[reaction_id]
        masks.append(mask)

    return masks, len(universe)


def pairwise_jaccard_from_masks(
    masks: list[int],
    progress_every: int = 100,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute only i < j and mirror into a full symmetric matrix."""
    n = len(masks)
    n_pairs = n * (n - 1) // 2

    matrix = np.zeros((n, n), dtype=np.float32)
    upper = np.empty(n_pairs, dtype=np.float32)

    k = 0
    started = time.time()
    for i in range(n - 1):
        left = masks[i]
        for j in range(i + 1, n):
            right = masks[j]
            union_count = (left | right).bit_count()
            if union_count == 0:
                distance = 0.0
            else:
                intersection_count = (left & right).bit_count()
                distance = 1.0 - (intersection_count / union_count)

            upper[k] = distance
            matrix[i, j] = distance
            matrix[j, i] = distance
            k += 1

        if progress_every > 0 and (
            (i + 1) % progress_every == 0 or i == n - 2
        ):
            elapsed = time.time() - started
            print(
                f"  Jaccard rows: {i + 1}/{n - 1}; "
                f"pairs={k}/{n_pairs}; elapsed={elapsed:.1f}s",
                flush=True,
            )

    return matrix, upper


def write_matrix_tsv_gz(
    path: Path,
    model_ids: list[str],
    matrix: np.ndarray,
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["model_id", *model_ids])
        for model_id, row in zip(model_ids, matrix):
            writer.writerow(
                [model_id, *[f"{float(value):.8f}" for value in row]]
            )


def compute_one_set(
    set_name: str,
    model_ids: list[str],
    mapping: dict[str, Path],
    output_root: Path,
    write_tsv: bool,
    progress_every: int,
) -> dict:
    started = time.time()
    set_dir = output_root / set_name
    set_dir.mkdir(parents=True, exist_ok=True)

    print(f"[{set_name}] Reading {len(model_ids)} SBML models", flush=True)

    reaction_sets: list[set[str]] = []
    reaction_counts: list[int] = []
    for idx, model_id in enumerate(model_ids, start=1):
        reactions = read_reaction_ids(mapping[model_id])
        reaction_sets.append(reactions)
        reaction_counts.append(len(reactions))
        if progress_every > 0 and (
            idx % progress_every == 0 or idx == len(model_ids)
        ):
            print(
                f"  parsed models: {idx}/{len(model_ids)}",
                flush=True,
            )

    masks, universe_size = reaction_sets_to_bitmasks(reaction_sets)
    # Release the large Python set objects before allocating/writing matrices.
    del reaction_sets

    print(
        f"[{set_name}] Reaction universe: {universe_size}; "
        f"computing upper triangle only",
        flush=True,
    )
    matrix, upper = pairwise_jaccard_from_masks(
        masks,
        progress_every=progress_every,
    )

    np.save(set_dir / "jaccard_reaction_distance.npy", matrix)
    np.save(set_dir / "jaccard_reaction_distance_upper.npy", upper)
    np.savez_compressed(
        set_dir / "jaccard_reaction_distance.npz",
        matrix=matrix,
        model_ids=np.asarray(model_ids, dtype=str),
    )

    with (set_dir / "models.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["index", "model_id", "n_reactions", "model_path"])
        for idx, model_id in enumerate(model_ids):
            writer.writerow([
                idx,
                model_id,
                reaction_counts[idx],
                str(mapping[model_id]),
            ])

    if write_tsv:
        print(f"[{set_name}] Writing compressed labeled TSV", flush=True)
        write_matrix_tsv_gz(
            set_dir / "jaccard_reaction_distance.tsv.gz",
            model_ids,
            matrix,
        )

    summary = {
        "set_name": set_name,
        "n_models": len(model_ids),
        "n_unique_pairs_computed": int(len(upper)),
        "matrix_shape": [int(matrix.shape[0]), int(matrix.shape[1])],
        "distance_definition": "1 - |reaction_ids_i intersection reaction_ids_j| / |reaction_ids_i union reaction_ids_j|",
        "reaction_universe_size": int(universe_size),
        "reaction_count_min": int(min(reaction_counts)),
        "reaction_count_max": int(max(reaction_counts)),
        "reaction_count_mean": float(np.mean(reaction_counts)),
        "jaccard_distance_min_off_diagonal": float(np.min(upper)),
        "jaccard_distance_max_off_diagonal": float(np.max(upper)),
        "jaccard_distance_mean_off_diagonal": float(np.mean(upper)),
        "jaccard_distance_median_off_diagonal": float(np.median(upper)),
        "elapsed_seconds": float(time.time() - started),
        "outputs": {
            "matrix_npy": str(set_dir / "jaccard_reaction_distance.npy"),
            "upper_triangle_npy": str(set_dir / "jaccard_reaction_distance_upper.npy"),
            "matrix_npz": str(set_dir / "jaccard_reaction_distance.npz"),
            "matrix_tsv_gz": (
                str(set_dir / "jaccard_reaction_distance.tsv.gz")
                if write_tsv
                else None
            ),
            "models_tsv": str(set_dir / "models.tsv"),
        },
    }

    with (set_dir / "summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(
        f"[{set_name}] done in {summary['elapsed_seconds']:.1f}s; "
        f"mean distance={summary['jaccard_distance_mean_off_diagonal']:.6f}",
        flush=True,
    )
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compute three 1375x1375 reaction-content Jaccard distance matrices "
            "for TARA Chile gapseq models."
        )
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=Path("~/tara_chile_metabolic_models"),
        help="Root of the tara_chile_metabolic_models checkout.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        required=True,
        help="Output directory for the three matrices.",
    )
    parser.add_argument(
        "--expected-models",
        type=int,
        default=1375,
        help="Required number of models in each set (default: 1375).",
    )
    parser.add_argument(
        "--no-tsv",
        action="store_true",
        help="Skip the labeled .tsv.gz matrices; .npy/.npz are always written.",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=100,
        help="Progress reporting interval in models/rows.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate directory counts and matching model IDs, then exit.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data_root = args.data_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()

    print("============================================================", flush=True)
    print("TARA Chile model Jaccard distance analysis", flush=True)
    print(f"Data root      : {data_root}", flush=True)
    print(f"Output root    : {output_root}", flush=True)
    print(f"Expected models: {args.expected_models}", flush=True)
    print("Distance basis : exact SBML reaction IDs", flush=True)
    print("Computation    : strict upper triangle only; mirrored afterward", flush=True)
    print("============================================================", flush=True)

    model_ids, mappings = validate_model_sets(
        data_root=data_root,
        expected_models=args.expected_models,
    )

    print(
        f"Validated all three sets: {len(model_ids)} matching model IDs "
        f"({model_ids[0]} ... {model_ids[-1]})",
        flush=True,
    )

    if args.validate_only:
        return 0

    output_root.mkdir(parents=True, exist_ok=True)

    with (output_root / "model_ids.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["index", "model_id"])
        for idx, model_id in enumerate(model_ids):
            writer.writerow([idx, model_id])

    summaries = []
    for set_name, _relative_path in DEFAULT_SETS:
        summaries.append(
            compute_one_set(
                set_name=set_name,
                model_ids=model_ids,
                mapping=mappings[set_name],
                output_root=output_root,
                write_tsv=not args.no_tsv,
                progress_every=args.progress_every,
            )
        )

    run_summary = {
        "n_models": len(model_ids),
        "n_unique_pairs_per_set": len(model_ids) * (len(model_ids) - 1) // 2,
        "sets": summaries,
    }
    with (output_root / "run_summary.json").open("w") as handle:
        json.dump(run_summary, handle, indent=2)
        handle.write("\n")

    print("============================================================", flush=True)
    print("All three Jaccard matrices completed successfully.", flush=True)
    print(f"Results: {output_root}", flush=True)
    print("============================================================", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
