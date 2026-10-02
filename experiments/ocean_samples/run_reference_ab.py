"""Stage 1: full-community reference scan plus product Stages A and B."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from time import perf_counter

from reframed.solvers.solution import Status

from common import yaml_dump_atomic, yaml_load
from misosoup.library.product_filter import (
    filter_exchange_candidates,
    filter_product_candidates,
)
from misosoup.library.product_reference import (
    constrain_full_community_lp,
    find_producible_exchanges,
    get_community_exchanges,
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

        all_exchange_reactions = sorted(
            get_community_exchanges(community)
        )
        scan_candidates, prefilter_audit = filter_exchange_candidates(
            community=community,
            exchange_reactions=all_exchange_reactions,
            keep_oxygen=args.keep_oxygen,
        )
        prefilter_excluded = {
            rid: info
            for rid, info in prefilter_audit.items()
            if not info["keep"]
        }

        print(
            "Product prefilter: "
            f"{len(all_exchange_reactions)} total exchanges -> "
            f"{len(scan_candidates)} selected for optimization; "
            f"{len(prefilter_excluded)} excluded before LP solves.",
            flush=True,
        )
        if not scan_candidates:
            raise RuntimeError(
                "Product prefilter removed every global exchange; "
                "there are no candidate products to optimize."
            )

        scan_checkpoint_path = stage_dir / "product_scan.yaml"
        scan_progress_path = stage_dir / "product_scan_partial.yaml"

        scan_configuration = {
            "checkpoint_version": 2,
            "minimal_growth": float(args.minimal_growth),
            "production_tolerance": float(
                args.production_tolerance
            ),
            "keep_oxygen": bool(args.keep_oxygen),
            "prefilter_before_optimization": True,
            "medium_file": str(
                args.medium_file.expanduser().resolve()
            ),
        }

        def scan_state_matches(payload: dict) -> bool:
            configuration = payload.get("configuration", {})
            return (
                int(configuration.get("checkpoint_version", -1)) == 2
                and float(configuration.get("minimal_growth", -1))
                == float(args.minimal_growth)
                and float(
                    configuration.get("production_tolerance", -1)
                )
                == float(args.production_tolerance)
                and bool(configuration.get("keep_oxygen", False))
                == bool(args.keep_oxygen)
                and bool(
                    configuration.get(
                        "prefilter_before_optimization",
                        False,
                    )
                )
                is True
                and str(configuration.get("medium_file", ""))
                == str(args.medium_file.expanduser().resolve())
                and list(
                    payload.get("candidate_exchange_reactions", [])
                )
                == list(scan_candidates)
            )

        if scan_checkpoint_path.exists():
            cached_scan = yaml_load(scan_checkpoint_path)
            cache_matches = (
                cached_scan.get("status") == "complete"
                and scan_state_matches(cached_scan)
            )
        else:
            cache_matches = False

        if cache_matches:
            unfiltered = {
                rid: float(value)
                for rid, value in cached_scan[
                    "unfiltered_max_secretion"
                ].items()
            }
            filtered = {
                rid: float(value)
                for rid, value in cached_scan[
                    "filtered_max_secretion"
                ].items()
            }
            filter_audit = cached_scan["filter_audit"]
            prefilter_audit = cached_scan["prefilter_audit"]
            scan_seconds = float(cached_scan.get("scan_seconds", 0.0))
            print(
                f"Product scan checkpoint reused: "
                f"{len(unfiltered)} producible -> "
                f"{len(filtered)} filtered candidates.",
                flush=True,
            )
        else:
            exchange_results = {}

            if scan_progress_path.exists():
                partial_scan = yaml_load(scan_progress_path)

                if scan_state_matches(partial_scan):
                    exchange_results = dict(
                        partial_scan.get("exchange_results", {})
                    )
                    restored_optimal = sum(
                        1
                        for record in exchange_results.values()
                        if record.get("status") == "optimal"
                        and record.get("maximum") is not None
                    )
                    print(
                        "Product scan partial checkpoint restored: "
                        f"{restored_optimal}/{len(scan_candidates)} "
                        "optimal exchanges already complete; "
                        f"{len(scan_candidates) - restored_optimal} "
                        "remaining.",
                        flush=True,
                    )
                else:
                    print(
                        "Existing product_scan_partial.yaml does not "
                        "match the current scan configuration/candidate "
                        "set and will be replaced.",
                        flush=True,
                    )

            def persist_scan_result(
                reaction_id: str,
                record: dict,
            ) -> None:
                exchange_results[reaction_id] = record

                completed_optimal = sum(
                    1
                    for value in exchange_results.values()
                    if value.get("status") == "optimal"
                    and value.get("maximum") is not None
                )
                attempted = len(exchange_results)
                solve_seconds = sum(
                    float(value.get("solve_seconds", 0.0))
                    for value in exchange_results.values()
                )

                yaml_dump_atomic(
                    {
                        "status": "partial",
                        "configuration": scan_configuration,
                        "candidate_exchange_reactions": list(
                            scan_candidates
                        ),
                        "number_total_exchange_reactions": len(
                            all_exchange_reactions
                        ),
                        "number_prefilter_candidates": len(
                            scan_candidates
                        ),
                        "number_prefilter_excluded": len(
                            prefilter_excluded
                        ),
                        "number_attempted": attempted,
                        "number_completed_optimal": completed_optimal,
                        "number_remaining_optimal": (
                            len(scan_candidates) - completed_optimal
                        ),
                        "cumulative_completed_solve_seconds": (
                            solve_seconds
                        ),
                        "exchange_results": exchange_results,
                    },
                    scan_progress_path,
                )

            t = perf_counter()
            unfiltered = find_producible_exchanges(
                community,
                tolerance=args.production_tolerance,
                exchange_reactions=scan_candidates,
                completed_results=exchange_results,
                progress_callback=persist_scan_result,
            )
            current_scan_wall_seconds = perf_counter() - t

            scan_seconds = sum(
                float(value.get("solve_seconds", 0.0))
                for value in exchange_results.values()
            )

            nonoptimal_results = {
                rid: value
                for rid, value in exchange_results.items()
                if value.get("status") != "optimal"
                or value.get("maximum") is None
            }
            if nonoptimal_results:
                yaml_dump_atomic(
                    {
                        "status": "partial",
                        "configuration": scan_configuration,
                        "candidate_exchange_reactions": list(
                            scan_candidates
                        ),
                        "number_total_exchange_reactions": len(
                            all_exchange_reactions
                        ),
                        "number_prefilter_candidates": len(
                            scan_candidates
                        ),
                        "number_prefilter_excluded": len(
                            prefilter_excluded
                        ),
                        "number_attempted": len(exchange_results),
                        "number_completed_optimal": (
                            len(exchange_results)
                            - len(nonoptimal_results)
                        ),
                        "number_remaining_optimal": len(
                            nonoptimal_results
                        ),
                        "cumulative_completed_solve_seconds": (
                            scan_seconds
                        ),
                        "exchange_results": exchange_results,
                    },
                    scan_progress_path,
                )
                raise RuntimeError(
                    "Product scan has non-optimal exchange solves; "
                    "partial checkpoint preserved so they can be "
                    "retried on the next run: "
                    + ", ".join(
                        sorted(nonoptimal_results)[:20]
                    )
                )

            filtered, filter_audit = filter_product_candidates(
                community=community,
                max_secretion=unfiltered,
                keep_oxygen=args.keep_oxygen,
            )

            final_scan_payload = {
                "status": "complete",
                "configuration": scan_configuration,
                "candidate_exchange_reactions": list(
                    scan_candidates
                ),
                "number_total_exchange_reactions": len(
                    all_exchange_reactions
                ),
                "number_prefilter_candidates": len(scan_candidates),
                "number_prefilter_excluded": len(prefilter_excluded),
                "number_attempted": len(exchange_results),
                "number_completed_optimal": sum(
                    1
                    for value in exchange_results.values()
                    if value.get("status") == "optimal"
                    and value.get("maximum") is not None
                ),
                "prefilter_audit": prefilter_audit,
                "exchange_results": exchange_results,
                "unfiltered_max_secretion": {
                    rid: float(value)
                    for rid, value in unfiltered.items()
                },
                "filtered_max_secretion": {
                    rid: float(value)
                    for rid, value in filtered.items()
                },
                "filter_audit": filter_audit,
                "scan_seconds": float(scan_seconds),
                "current_run_scan_wall_seconds": float(
                    current_scan_wall_seconds
                ),
            }

            yaml_dump_atomic(
                final_scan_payload,
                scan_checkpoint_path,
            )
            yaml_dump_atomic(
                {
                    **final_scan_payload,
                    "status": "complete",
                },
                scan_progress_path,
            )

            print(
                f"Product scan complete: {len(unfiltered)} producible -> "
                f"{len(filtered)} filtered candidates; "
                f"cumulative completed-solve time="
                f"{scan_seconds:.3f}s; final checkpoint written.",
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
                "number_total_exchange_reactions": len(
                    all_exchange_reactions
                ),
                "number_prefilter_candidates": len(scan_candidates),
                "number_prefilter_excluded": len(prefilter_excluded),
                "number_unfiltered_producible_exchanges": len(unfiltered),
                "number_filtered_products": len(filtered),
                "unfiltered_max_secretion": {
                    rid: float(value) for rid, value in unfiltered.items()
                },
                "filtered_max_secretion": {
                    rid: float(value) for rid, value in filtered.items()
                },
                "prefilter_audit": prefilter_audit,
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
