"""Run the validated product-preservation workflow for one ocean sample."""

from __future__ import annotations

import argparse
import copy
import json
import os
import traceback
from pathlib import Path
from time import perf_counter

from reframed.solvers.solution import Status
from reframed.solvers.solver import Parameter

from common import (
    read_medium_spec,
    resolve_medium_for_community,
    yaml_dump_atomic,
    yaml_load,
)
from misosoup.library.product_filter import filter_product_candidates
from misosoup.library.product_reference import (
    constrain_full_community,
    find_producible_exchanges,
)
from misosoup.library.minimal_product_communities import (
    findMinimalProductCommunities,
)
from misosoup.library.product_selection import (
    getMaxProduct,
    getSelectedProductsFromProductMaximization,
    maximizeProducts,
    solvepFBAUsingFixProducts,
)
from misosoup.library.readwrite import load_models
from misosoup.reframed.layered_community import LayeredCommunity


def elapsed(start: float) -> float:
    return perf_counter() - start


def solver_params(prod_tol: float, int_tol: float) -> dict:
    return {
        Parameter.OPTIMALITY_TOL: prod_tol,
        Parameter.FEASIBILITY_TOL: prod_tol,
        Parameter.INT_FEASIBILITY_TOL: int_tol,
    }


def write_json_atomic(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
    tmp.replace(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-index", required=True, type=int)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=0.20)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--product-retention", type=float, default=0.90)
    parser.add_argument("--lexicographic-tolerance", type=float, default=1e-5)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument(
        "--keep-oxygen",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument("--max-minimal-communities", type=int, default=100)
    return parser.parse_args()


def map_stage_d(stage_d: dict, solver_to_mag: dict[str, str]) -> dict:
    converted = copy.deepcopy(stage_d)
    for community in converted["communities"]:
        solver_ids = list(community["organisms"])
        community["solver_organisms"] = solver_ids
        community["organisms"] = [solver_to_mag[x] for x in solver_ids]
        community["mags"] = list(community["organisms"])
        solver_growth = dict(community.get("organism_growth", {}))
        community["solver_organism_growth"] = solver_growth
        community["organism_growth"] = {
            solver_to_mag[x]: float(v)
            for x, v in solver_growth.items()
        }
    return converted


def main() -> None:
    args = parse_args()

    if not 0 < args.alpha <= 1:
        raise ValueError("--alpha must be >0 and <=1")
    if args.minimal_growth <= 0:
        raise ValueError("--minimal-growth must be >0")
    if not 0 < args.product_retention <= 1:
        raise ValueError("--product-retention must be >0 and <=1")
    if args.lexicographic_tolerance < 0 or args.production_tolerance < 0:
        raise ValueError("Tolerances must be non-negative")
    if args.uptake_bound >= 0:
        raise ValueError("--uptake-bound must be negative")
    if args.max_minimal_communities < 1:
        raise ValueError("--max-minimal-communities must be >=1")

    manifest = yaml_load(args.manifest.expanduser().resolve())
    samples = manifest["samples"]
    if not 0 <= args.sample_index < len(samples):
        raise IndexError(
            f"sample index {args.sample_index} outside 0..{len(samples)-1}"
        )

    sample = samples[args.sample_index]
    sample_id = str(sample["sample_id"])
    output_root = args.output_root.expanduser().resolve()
    sample_dir = output_root / "samples" / f"{args.sample_index:03d}"
    sample_dir.mkdir(parents=True, exist_ok=True)
    medium_file = args.medium_file.expanduser().resolve()
    total_start = perf_counter()

    write_json_atomic(
        {
            "sample_index": args.sample_index,
            "sample_id": sample_id,
            "status": "running",
            "pid": os.getpid(),
        },
        sample_dir / "status.json",
    )

    try:
        print("=" * 78, flush=True)
        print(
            f"Ocean sample {args.sample_index}: {sample_id}; "
            f"MAGs={sample['number_modeled_mags']}",
            flush=True,
        )
        print("=" * 78, flush=True)

        model_paths = [
            Path(entry["model_path"]).expanduser().resolve()
            for entry in sample["mags"]
        ]
        missing = [str(p) for p in model_paths if not p.exists()]
        if missing:
            raise FileNotFoundError(
                "Missing model files: " + ", ".join(missing[:20])
            )

        t = perf_counter()
        models = load_models([str(p) for p in model_paths])
        model_loading = elapsed(t)

        if len(models) != len(sample["mags"]):
            raise RuntimeError("Loaded-model count differs from manifest")

        solver_to_mag = {}
        mag_to_solver = {}
        for position, (model, entry) in enumerate(
            zip(models, sample["mags"])
        ):
            solver_id = f"org_{position:05d}"
            mag_id = str(entry["mag_id"])
            model.id = solver_id
            solver_to_mag[solver_id] = mag_id
            mag_to_solver[mag_id] = solver_id

        params = solver_params(
            args.production_tolerance,
            args.integer_tolerance,
        )

        t = perf_counter()
        community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_reference",
            models,
            copy_models=False,
            params=params,
        )
        community_build = elapsed(t)

        constrain_full_community(
            community,
            minimal_growth=args.minimal_growth,
        )

        medium_entries, medium_file_audit = read_medium_spec(
            medium_file,
            default_uptake_bound=args.uptake_bound,
        )
        medium, medium_audit = resolve_medium_for_community(
            community,
            medium_entries,
            args.uptake_bound,
        )
        community.setup_medium(medium)

        t = perf_counter()
        feasibility = community.check_feasibility(["community_growth"])
        feasibility_time = elapsed(t)
        if feasibility.status != Status.OPTIMAL:
            raise RuntimeError(
                "Full sample community is not feasible under the selected medium "
                f"at minimal_growth={args.minimal_growth}; "
                f"solver status={feasibility.status}"
            )

        print(
            f"Medium matched {medium_audit['number_matched_tokens']} / "
            f"{medium_audit['number_tokens']} tokens",
            flush=True,
        )

        t = perf_counter()
        unfiltered = find_producible_exchanges(
            community,
            tolerance=args.production_tolerance,
        )
        scan_time = elapsed(t)

        filtered, filter_audit = filter_product_candidates(
            community=community,
            max_secretion=unfiltered,
            keep_oxygen=args.keep_oxygen,
        )
        excluded = {
            rid: info
            for rid, info in filter_audit.items()
            if not info["keep"]
        }
        missing_formula = {
            rid: info
            for rid, info in filter_audit.items()
            if info["reason"] == "missing_formula_kept_for_review"
        }

        print(
            f"Product scan: {len(unfiltered)} producible -> "
            f"{len(filtered)} filtered candidates in {scan_time:.3f}s",
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
        stage_a_time = elapsed(t)

        t = perf_counter()
        stage_b = maximizeProducts(
            community=community,
            stage_a=stage_a,
        )
        stage_b_time = elapsed(t)

        t = perf_counter()
        product_selection = getSelectedProductsFromProductMaximization(
            stage_b=stage_b
        )
        stage_b2_time = elapsed(t)

        t = perf_counter()
        stage_c_community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_stage_c",
            models,
            copy_models=False,
            params=params,
        )
        stage_c = solvepFBAUsingFixProducts(
            community=stage_c_community,
            stage_b=stage_b,
            product_selection=product_selection,
            medium=medium,
            minimal_growth=args.minimal_growth,
            lexicographic_tolerance=args.lexicographic_tolerance,
            tolerance=args.production_tolerance,
            check_feasibility=True,
        )
        stage_c_time = elapsed(t)

        if set(stage_c["selected_products"]) != set(
            product_selection["selected_products"]
        ):
            raise RuntimeError("Stage C product set differs from Stage B")

        t = perf_counter()
        stage_d_community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_stage_d",
            models,
            copy_models=False,
            params=params,
        )
        stage_d_raw = findMinimalProductCommunities(
            community=stage_d_community,
            medium=medium,
            reference_products=stage_c["reference_products"],
            minimal_growth=args.minimal_growth,
            product_retention=args.product_retention,
            tolerance=args.production_tolerance,
            max_communities=args.max_minimal_communities,
        )
        stage_d_time = elapsed(t)
        stage_d = map_stage_d(stage_d_raw, solver_to_mag)

        total_time = elapsed(total_start)

        selected_info = {
            rid: {key: float(value) for key, value in info.items()}
            for rid, info in product_selection[
                "selected_product_info"
            ].items()
        }
        reference_products = {
            rid: {key: float(value) for key, value in info.items()}
            for rid, info in stage_c["reference_products"].items()
        }

        result = {
            "status": "ok",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": sample_id,
                "number_present_mags_in_matrix": int(
                    sample["number_present_mags_in_matrix"]
                ),
                "number_modeled_mags": int(
                    sample["number_modeled_mags"]
                ),
                "members": [
                    {
                        "mag_id": str(entry["mag_id"]),
                        "abundance": float(entry["abundance"]),
                        "model_path": str(entry["model_path"]),
                        "solver_id": mag_to_solver[str(entry["mag_id"])],
                    }
                    for entry in sample["mags"]
                ],
            },
            "configuration": {
                "alpha": float(args.alpha),
                "minimal_growth": float(args.minimal_growth),
                "product_retention": float(args.product_retention),
                "lexicographic_tolerance": float(
                    args.lexicographic_tolerance
                ),
                "production_tolerance": float(args.production_tolerance),
                "integer_tolerance": float(args.integer_tolerance),
                "keep_oxygen": bool(args.keep_oxygen),
                "uptake_bound": float(args.uptake_bound),
                "max_minimal_communities": int(
                    args.max_minimal_communities
                ),
                "abundance_used_as_metabolic_weight": False,
            },
            "medium": {
                "file": str(medium_file),
                **medium_file_audit,
                **medium_audit,
                "resolved_exchange_bounds": {
                    rid: float(bound)
                    for rid, bound in medium.items()
                },
            },
            "individual_product_maxima": {
                "number_unfiltered_producible_exchanges": len(unfiltered),
                "number_filtered_products": len(filtered),
                "unfiltered_max_secretion": {
                    rid: float(v) for rid, v in unfiltered.items()
                },
                "filtered_max_secretion": {
                    rid: float(v) for rid, v in filtered.items()
                },
                "excluded_products": excluded,
                "missing_formula_products": missing_formula,
            },
            "stage_a": {
                "k_star": int(stage_a["k_star"]),
                "selected_products": list(
                    stage_a["stage_a_selected_products"]
                ),
            },
            "stage_b": {
                "q_star": float(stage_b["q_star"]),
                "q_star_from_values": float(
                    stage_b["q_star_from_values"]
                ),
                "q_star_difference": float(
                    stage_b["q_star_difference"]
                ),
                "selected_products": list(
                    product_selection["selected_products"]
                ),
                "selected_product_info": selected_info,
            },
            "stage_c": {
                "selected_products": list(stage_c["selected_products"]),
                "stage_b_objective_floor": float(
                    stage_c["stage_b_objective_floor"]
                ),
                "stage_c_q_sum": float(stage_c["stage_c_q_sum"]),
                "lexicographic_tolerance": float(
                    stage_c["lexicographic_tolerance"]
                ),
                "pfba_objective": float(stage_c["pfba_objective"]),
                "reference_products": reference_products,
            },
            "stage_d": stage_d,
            "timing_seconds": {
                "model_loading": model_loading,
                "community_build": community_build,
                "full_community_feasibility": feasibility_time,
                "individual_product_scan": scan_time,
                "stage_a": stage_a_time,
                "stage_b1": stage_b_time,
                "stage_b2": stage_b2_time,
                "stage_c": stage_c_time,
                "stage_d": stage_d_time,
                "total": total_time,
            },
        }

        yaml_dump_atomic(result, sample_dir / "result.yaml")
        write_json_atomic(
            {
                "sample_index": int(args.sample_index),
                "sample_id": sample_id,
                "status": "ok",
                "number_modeled_mags": int(
                    sample["number_modeled_mags"]
                ),
                "k_star": int(stage_a["k_star"]),
                "q_star": float(stage_b["q_star"]),
                "number_selected_products": len(
                    stage_c["selected_products"]
                ),
                "pfba_objective": float(stage_c["pfba_objective"]),
                "minimum_size": int(stage_d["minimum_size"]),
                "number_minimum_communities": int(
                    stage_d["number_communities"]
                ),
                "enumeration_complete": bool(
                    stage_d["enumeration_complete"]
                ),
                "total_seconds": float(total_time),
            },
            sample_dir / "status.json",
        )

        print(
            f"COMPLETED sample {args.sample_index}: "
            f"K*={stage_a['k_star']} "
            f"Q*={stage_b['q_star']:.12g} "
            f"N_min={stage_d['minimum_size']} "
            f"communities={stage_d['number_communities']} "
            f"total={total_time:.3f}s",
            flush=True,
        )

    except Exception as error:
        total_time = elapsed(total_start)
        trace = traceback.format_exc()
        failure = {
            "status": "failed",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": sample_id,
                "number_present_mags_in_matrix": int(
                    sample["number_present_mags_in_matrix"]
                ),
                "number_modeled_mags": int(
                    sample["number_modeled_mags"]
                ),
            },
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": trace,
            "timing_seconds": {"total": total_time},
        }
        yaml_dump_atomic(failure, sample_dir / "failure.yaml")
        write_json_atomic(
            {
                "sample_index": int(args.sample_index),
                "sample_id": sample_id,
                "status": "failed",
                "error_type": type(error).__name__,
                "error": str(error),
                "total_seconds": float(total_time),
            },
            sample_dir / "status.json",
        )
        with (sample_dir / "traceback.txt").open(
            "w", encoding="utf-8"
        ) as handle:
            handle.write(trace)
        print(trace, flush=True)
        raise


if __name__ == "__main__":
    main()
