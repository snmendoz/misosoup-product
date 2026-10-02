"""Stage 2: reconstruct the reference community and solve pure-LP pFBA."""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from common import yaml_dump_atomic, yaml_load
from misosoup.library.product_selection import solvepFBAUsingFixProducts
from misosoup.reframed.layered_community import LayeredCommunity
from staged_common import (
    get_sample,
    load_sample_models,
    mark_complete,
    resolve_medium,
    sample_dir,
    solver_params,
    write_stage_failure,
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-index", required=True, type=int)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--lexicographic-tolerance", type=float, default=1e-5)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    return parser.parse_args()


def main():
    args = parse_args()
    _, sample = get_sample(args.manifest, args.sample_index)
    root = sample_dir(args.output_root, args.sample_index)
    stage_dir = root / "02_stage_c"
    stage_dir.mkdir(parents=True, exist_ok=True)
    total_start = perf_counter()

    try:
        ab_path = root / "01_reference_ab" / "reference_ab.yaml"
        if not ab_path.exists():
            raise FileNotFoundError(f"Missing Stage A/B checkpoint: {ab_path}")
        ab = yaml_load(ab_path)
        if ab.get("status") != "complete":
            raise RuntimeError("Stage A/B checkpoint is not complete.")

        models, _, _, load_seconds = load_sample_models(sample)
        params = solver_params(args.production_tolerance, args.integer_tolerance)

        t = perf_counter()
        community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_stage_c",
            models,
            copy_models=False,
            params=params,
        )
        build_seconds = perf_counter() - t

        medium, medium_audit, medium_file_audit = resolve_medium(
            community,
            args.medium_file,
            args.uptake_bound,
        )

        t = perf_counter()
        stage_c = solvepFBAUsingFixProducts(
            community=community,
            stage_b=ab["stage_b"],
            product_selection=ab["product_selection"],
            medium=medium,
            minimal_growth=args.minimal_growth,
            lexicographic_tolerance=args.lexicographic_tolerance,
            tolerance=args.production_tolerance,
            check_feasibility=True,
        )
        stage_c_seconds = perf_counter() - t

        checkpoint = {
            "status": "complete",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": str(sample["sample_id"]),
                "number_modeled_mags": int(sample["number_modeled_mags"]),
            },
            "configuration": {
                "minimal_growth": float(args.minimal_growth),
                "lexicographic_tolerance": float(
                    args.lexicographic_tolerance
                ),
                "production_tolerance": float(args.production_tolerance),
                "integer_tolerance": float(args.integer_tolerance),
                "uptake_bound": float(args.uptake_bound),
                "medium_file": str(args.medium_file.expanduser().resolve()),
            },
            "medium": {
                **medium_file_audit,
                **medium_audit,
            },
            "stage_c": stage_c,
            "timing_seconds": {
                "model_loading": float(load_seconds),
                "community_build": float(build_seconds),
                "stage_c": float(stage_c_seconds),
                "total": float(perf_counter() - total_start),
            },
        }

        yaml_dump_atomic(checkpoint, stage_dir / "stage_c.yaml")
        mark_complete(stage_dir)
        print(
            f"STAGE C COMPLETE in {perf_counter() - total_start:.3f}s.",
            flush=True,
        )

    except Exception as error:
        write_stage_failure(
            stage_name="stage_c",
            stage_dir=stage_dir,
            sample_index=args.sample_index,
            sample_id=str(sample["sample_id"]),
            error=error,
            elapsed_seconds=perf_counter() - total_start,
            sample_root=root,
        )
        raise


if __name__ == "__main__":
    main()
