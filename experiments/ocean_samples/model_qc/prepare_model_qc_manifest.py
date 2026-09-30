"""Build a strict MAG-to-model manifest for individual model QC."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import (
    build_model_index,
    discover_model_files,
    read_abundance_matrix,
    resolve_mag_model,
    yaml_dump_atomic,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True, type=Path)
    parser.add_argument("--models-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--expected-samples", type=int, default=159)
    parser.add_argument("--expected-mags", type=int, default=1375)
    parser.add_argument(
        "--mag-ids",
        default=None,
        help=(
            "Optional comma-separated MAG IDs. "
            "If omitted, all MAGs in the abundance matrix are checked."
        ),
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    matrix = read_abundance_matrix(
        args.matrix,
        expected_samples=args.expected_samples,
        expected_mags=args.expected_mags,
    )

    all_mag_ids = [str(value) for value in matrix.columns]

    if args.mag_ids:
        requested = [
            value.strip()
            for value in args.mag_ids.split(",")
            if value.strip()
        ]

        missing = sorted(set(requested) - set(all_mag_ids))

        if missing:
            raise ValueError(
                "Requested MAG IDs are not present in the abundance matrix: "
                + ", ".join(missing)
            )

        mag_ids = requested
    else:
        mag_ids = all_mag_ids

    model_paths = discover_model_files(args.models_dir)
    exact_index, normalized_index = build_model_index(model_paths)

    models = []
    used_paths = {}

    for index, mag_id in enumerate(mag_ids):
        model_path = resolve_mag_model(
            mag_id=mag_id,
            exact_index=exact_index,
            normalized_index=normalized_index,
        )

        model_key = str(model_path)

        if model_key in used_paths:
            raise ValueError(
                "Two MAG IDs map to the same model file: "
                f"{used_paths[model_key]} and {mag_id} -> {model_path}"
            )

        used_paths[model_key] = mag_id

        models.append(
            {
                "index": int(index),
                "mag_id": mag_id,
                "model_path": model_key,
            }
        )

    manifest = {
        "manifest_version": 1,
        "source": {
            "abundance_matrix": str(args.matrix.expanduser().resolve()),
            "models_dir": str(args.models_dir.expanduser().resolve()),
        },
        "number_models": len(models),
        "number_models_discovered": len(model_paths),
        "subset_requested": bool(args.mag_ids),
        "models": models,
    }

    yaml_dump_atomic(
        manifest,
        args.output.expanduser().resolve(),
    )

    print("=" * 72)
    print("Individual model QC manifest created")
    print("=" * 72)
    print(f"Models discovered : {len(model_paths)}")
    print(f"Models selected   : {len(models)}")
    print(f"Subset requested  : {bool(args.mag_ids)}")
    print(f"Output            : {args.output.expanduser().resolve()}")
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
