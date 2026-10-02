"""Stage 1: full-community reference scan plus product Stages A and B."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from time import perf_counter

from reframed.solvers.solution import Status

from common import yaml_dump_atomic
from misosoup.library.product_filter import filter_product_candidates
from misosoup.library.product_reference import (
    constrain_full_community_lp,
    find_producible_exchanges,
)
from misosoup.library.product_selection import (
    getMaxProduct,
    getSelectedProductsFromProductMaximization,
    maximizeProducts,
)
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
    parser.add_argument("--alpha", type=float, default=0.20)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument(
        "--keep-oxygen",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    return parser.parse_args()


def write_product_csv(path: Path, values: dict) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["exchange_reaction", "maximum_secretion"])
        for rid, value in sorted(values.items()):
            writer.writerow([rid, float(value)])


def main():
    args = parse_args()
    _, sample = get_sample(args.manifest, args.sample_index)
    root = sample_dir(args.output_root, args.sample_index)
    stage_dir = root / "01_reference_ab"
    stage_dir.mkdir(parents=True, exist_ok=True)
    total_start = perf_counter()

    try:
        print("=" * 78, flush=True)
        print(
            f"REFERENCE + A/B: sample={sample['sample_id']} "
            f"MAGs={sample['number_modeled_mags']}",
            flush=True,
        )
        print("=" * 78, flush=True)

        models, _, _, load_seconds = load_sample_models(sample)
        params = solver_params(args.production_tolerance, args.integer_tolerance)

        t = perf_counter()
        community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_reference_ab",
            models,
            copy_models=False,
            params=params,
        )
        build_seconds = perf_counter() - t

        constrain_full_community_lp(
            community,
            minimal_growth=args.minimal_growth,
        )
        medium, medium_audit, medium_file_audit = resolve_medium(
            community,
            args.medium_file,
            args.uptake_bound,
        )
        community.setup_medium(medium)

        print(
            f"Medium matched {medium_audit['number_matched_tokens']} / "
            f"{medium_audit['number_tokens']} tokens.",
            flush=True,
        )

        t = perf_counter()
        feasibility = community.check_feasibility(["community_growth"])
        feasibility_seconds = perf_counter() - t
        print(
            f"Reference LP feasibility: {feasibility.status} "
            f"in {feasibility_seconds:.3f}s.",
            flush=True,
        )
        if feasibility.status != Status.OPTIMAL:
            raise RuntimeError(
                f"Full reference community is not feasible: {feasibility.status}"
            )

        t = perf_counter()
        unfiltered = find_producible_exchanges(
            community,
            tolerance=args.production_tolerance,
        )
        scan_seconds = perf_counter() - t

        filtered, filter_audit = filter_product_candidates(
            community=community,
            max_secretion=unfiltered,
            keep_oxygen=args.keep_oxygen,
        )
        print(
            f"Product scan complete: {len(unfiltered)} producible -> "
            f"{len(filtered)} filtered candidates in {scan_seconds:.3f}s.",
            flush=True,
        )

        t = perf_counter()
        stage_a = getMaxProduct(
            community=community,
            max_secretion=filtered,
            fraction=args.alpha,
            tolerance=args.production_tolerance,
            integer_tolerance=args.integer_tolerance,
        )
        stage_a_seconds = perf_counter() - t

        t = perf_counter()
        stage_b = maximizeProducts(
            community=community,
            stage_a=stage_a,
        )
        stage_b_seconds = perf_counter() - t

        t = perf_counter()
        selection = getSelectedProductsFromProductMaximization(stage_b=stage_b)
        stage_b2_seconds = perf_counter() - t

        stage_b_checkpoint = {
            "products": {
                rid: float(value)
                for rid, value in stage_b["products"].items()
            },
            "thresholds": {
                rid: float(value)
                for rid, value in stage_b["thresholds"].items()
            },
            "fraction": float(stage_b["fraction"]),
            "tolerance": float(stage_b["tolerance"]),
            "integer_tolerance": float(stage_b["integer_tolerance"]),
            "k_star": int(stage_b["k_star"]),
            "q_star": float(stage_b["q_star"]),
            "q_star_from_values": float(stage_b["q_star_from_values"]),
            "q_star_difference": float(stage_b["q_star_difference"]),
        }

        checkpoint = {
            "status": "complete",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": str(sample["sample_id"]),
                "number_modeled_mags": int(sample["number_modeled_mags"]),
            },
            "configuration": {
                "alpha": float(args.alpha),
                "minimal_growth": float(args.minimal_growth),
                "production_tolerance": float(args.production_tolerance),
                "integer_tolerance": float(args.integer_tolerance),
                "uptake_bound": float(args.uptake_bound),
                "keep_oxygen": bool(args.keep_oxygen),
                "medium_file": str(args.medium_file.expanduser().resolve()),
            },
            "medium": {
                **medium_file_audit,
                **medium_audit,
            },
            "individual_product_maxima": {
                "number_unfiltered_producible_exchanges": len(unfiltered),
                "number_filtered_products": len(filtered),
                "unfiltered_max_secretion": {
                    rid: float(value) for rid, value in unfiltered.items()
                },
                "filtered_max_secretion": {
                    rid: float(value) for rid, value in filtered.items()
                },
                "filter_audit": filter_audit,
            },
            "stage_a": {
                "k_star": int(stage_a["k_star"]),
                "fraction": float(stage_a["fraction"]),
                "stage_a_selected_products": list(
                    stage_a["stage_a_selected_products"]
                ),
            },
            "stage_b": stage_b_checkpoint,
            "product_selection": {
                "selected_products": list(selection["selected_products"]),
                "selected_product_info": {
                    rid: {
                        key: float(value)
                        for key, value in info.items()
                    }
                    for rid, info in selection["selected_product_info"].items()
                },
            },
            "timing_seconds": {
                "model_loading": float(load_seconds),
                "community_build": float(build_seconds),
                "feasibility": float(feasibility_seconds),
                "individual_product_scan": float(scan_seconds),
                "stage_a": float(stage_a_seconds),
                "stage_b1": float(stage_b_seconds),
                "stage_b2": float(stage_b2_seconds),
                "total": float(perf_counter() - total_start),
            },
        }

        yaml_dump_atomic(checkpoint, stage_dir / "reference_ab.yaml")
        write_product_csv(stage_dir / "producible_exchanges.csv", unfiltered)
        write_product_csv(stage_dir / "filtered_products.csv", filtered)
        mark_complete(stage_dir)

        print(
            f"REFERENCE + A/B COMPLETE in "
            f"{perf_counter() - total_start:.3f}s.",
            flush=True,
        )

    except Exception as error:
        write_stage_failure(
            stage_name="reference_ab",
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
