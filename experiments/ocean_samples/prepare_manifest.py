
"""Prepare a validated sample-to-model manifest for the ocean experiment."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import (
    build_model_index,
    discover_model_files,
    read_abundance_matrix,
    read_explicit_model_map,
    resolve_mag_model,
    yaml_dump_atomic,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a 159-sample manifest from a sample x MAG abundance matrix.")
    parser.add_argument("--matrix", required=True, type=Path, help="159 x 1375 abundance matrix; first column must contain sample IDs.")
    parser.add_argument("--models-dir", required=True, type=Path, help="Directory containing gapseq SBML models.")
    parser.add_argument("--output", required=True, type=Path, help="Manifest YAML to create.")
    parser.add_argument("--presence-threshold", type=float, default=0.0, help="A MAG is present when abundance > threshold. Default: 0.0.")
    parser.add_argument("--expected-samples", type=int, default=159, help="Expected number of samples. Default: 159.")
    parser.add_argument("--expected-mags", type=int, default=1375, help="Expected number of MAG columns. Default: 1375.")
    parser.add_argument("--model-map", type=Path, default=None, help="Optional explicit MAG-to-model table with columns mag_id and model_path.")
    parser.add_argument("--allow-missing-models", action="store_true", help="Skip present MAGs without a model instead of failing. Useful for temporary/pilot model sets; not recommended for the final analysis.")
    parser.add_argument("--max-mags-per-sample", type=int, default=None, help="Optional cap on modeled MAGs per sample. When set, the highest-abundance resolvable MAGs are selected first. Use 4 for the temporary pilot.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matrix_path = args.matrix.expanduser().resolve()
    models_dir = args.models_dir.expanduser().resolve()
    output_path = args.output.expanduser().resolve()

    if args.max_mags_per_sample is not None and args.max_mags_per_sample < 1:
        raise ValueError("--max-mags-per-sample must be >= 1 or omitted.")

    print(f"Reading abundance matrix: {matrix_path}", flush=True)
    frame = read_abundance_matrix(
        matrix_path,
        expected_samples=args.expected_samples,
        expected_mags=args.expected_mags,
    )
    print(
        f"Matrix shape after orientation check: "
        f"{frame.shape[0]} samples x {frame.shape[1]} MAGs",
        flush=True,
    )
    print(
        "Detected matrix orientation: "
        f"{frame.attrs.get('detected_orientation', 'unknown')}",
        flush=True,
    )
    if frame.attrs.get("mag_id_column"):
        print(
            "Detected MAG ID column: "
            f"{frame.attrs['mag_id_column']}",
            flush=True,
        )
    if frame.attrs.get("sample_id_column"):
        print(
            "Detected sample ID column: "
            f"{frame.attrs['sample_id_column']}",
            flush=True,
        )

    model_paths = discover_model_files(models_dir)
    exact_index, normalized_index = build_model_index(model_paths)
    explicit_map = read_explicit_model_map(args.model_map)
    print(f"Discovered {len(model_paths)} SBML models below {models_dir}", flush=True)

    samples = []
    unresolved_present_mags: dict[str, list[str]] = {}
    total_present_memberships = 0
    unique_present_mags: set[str] = set()
    resolved_models: set[str] = set()
    model_to_mag: dict[str, str] = {}

    for sample_index, (sample_id, row) in enumerate(frame.iterrows()):
        present = row[row > args.presence_threshold].sort_values(
            ascending=False
        )
        mag_entries = []
        unresolved = []
        for mag_id, abundance in present.items():
            mag_id = str(mag_id)
            unique_present_mags.add(mag_id)

            try:
                model_path = resolve_mag_model(
                    mag_id=mag_id,
                    exact_index=exact_index,
                    normalized_index=normalized_index,
                    explicit_map=explicit_map,
                )
            except FileNotFoundError:
                unresolved.append(mag_id)
                if args.allow_missing_models:
                    continue
                raise

            model_key = str(model_path)
            previous_mag = model_to_mag.get(model_key)

            if previous_mag is not None and previous_mag != mag_id:
                raise ValueError(
                    "Two abundance-matrix MAG IDs resolve to the same model: "
                    f"{previous_mag!r} and {mag_id!r} -> {model_path}"
                )

            model_to_mag[model_key] = mag_id
            resolved_models.add(model_key)
            mag_entries.append(
                {
                    "mag_id": mag_id,
                    "abundance": float(abundance),
                    "model_path": str(model_path),
                }
            )

            if (
                args.max_mags_per_sample is not None
                and len(mag_entries) >= args.max_mags_per_sample
            ):
                break

        if unresolved:
            unresolved_present_mags[str(sample_id)] = unresolved

        if not mag_entries:
            raise ValueError(
                f"Sample {sample_id!r} has no modeled MAGs above presence "
                f"threshold {args.presence_threshold} using the currently "
                "available model set."
            )

        total_present_memberships += len(mag_entries)

        samples.append(
            {
                "index": int(sample_index),
                "sample_id": str(sample_id),
                "number_present_mags_in_matrix": int(len(present)),
                "number_modeled_mags": int(len(mag_entries)),
                "mags": mag_entries,
            }
        )

    manifest = {
        "manifest_version": 1,
        "source": {
            "abundance_matrix": str(matrix_path),
            "models_dir": str(models_dir),
            "model_map": str(args.model_map.expanduser().resolve()) if args.model_map is not None else None,
        },
        "presence_rule": {
            "operator": ">",
            "threshold": float(args.presence_threshold),
            "abundance_used_as_metabolic_weight": False,
        },
        "matrix": {"number_samples": int(frame.shape[0]), "number_mag_columns": int(frame.shape[1])},
        "mapping": {
            "number_model_files_discovered": len(model_paths),
            "number_unique_present_mags": len(unique_present_mags),
            "number_unique_model_files_used": len(resolved_models),
            "number_sample_mag_memberships": total_present_memberships,
            "allow_missing_models": bool(args.allow_missing_models),
            "max_mags_per_sample": args.max_mags_per_sample,
            "selection_rule_when_capped": "highest_abundance_resolvable_MAGs_first",
            "unresolved_present_mags": unresolved_present_mags,
        },
        "samples": samples,
    }

    yaml_dump_atomic(manifest, output_path)

    print()
    print("=" * 72)
    print("Ocean-sample manifest created")
    print("=" * 72)
    print(f"Output                    : {output_path}")
    print(f"Samples                   : {len(samples)}")
    print(f"MAG columns               : {frame.shape[1]}")
    print(f"Unique present MAGs       : {len(unique_present_mags)}")
    print(f"Unique model files used   : {len(resolved_models)}")
    print(f"Sample-MAG memberships    : {total_present_memberships}")
    print(
        "Modeled MAGs per sample   : "
        f"min={min(s['number_modeled_mags'] for s in samples)}, "
        f"max={max(s['number_modeled_mags'] for s in samples)}, "
        f"mean={total_present_memberships / len(samples):.2f}"
    )
    print(f"Max MAGs per sample       : {args.max_mags_per_sample}")
    print(f"Missing mappings          : {sum(len(v) for v in unresolved_present_mags.values())}")
    print("=" * 72)
    print(flush=True)


if __name__ == "__main__":
    main()
