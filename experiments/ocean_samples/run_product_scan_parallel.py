"""Parallel Product Scan for one ocean-community sample.

Outer parallelism is supplied by a SLURM array across samples.  Inner
parallelism is supplied by forked SciPy/HiGHS workers across exchanges.
No Gurobi environment is created in this stage.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from time import perf_counter

from common import yaml_dump_atomic, yaml_load
from misosoup.library.highs_product_scan import (
    build_product_scan_lp,
    check_product_scan_feasibility,
    run_parallel_product_scan,
)
from misosoup.library.product_filter import (
    filter_exchange_candidates,
    filter_product_candidates,
)
from misosoup.library.product_reference import get_community_exchanges
from misosoup.reframed.layered_community import LayeredCommunity
from staged_common import (
    get_sample,
    load_sample_models,
    resolve_medium,
    sample_dir,
)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--sample-index", required=True, type=int)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument(
        "--workers",
        type=int,
        default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")),
    )
    parser.add_argument(
        "--keep-oxygen",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    return parser.parse_args()


def main():
    args = parse_args()
    _, sample = get_sample(args.manifest, args.sample_index)
    root = sample_dir(args.output_root, args.sample_index)
    stage_dir = root / "01_reference_ab"
    stage_dir.mkdir(parents=True, exist_ok=True)
    total_start = perf_counter()

    print("=" * 78, flush=True)
    print(
        f"PRODUCT SCAN: sample={sample['sample_id']} "
        f"MAGs={sample['number_modeled_mags']} "
        f"workers={args.workers}",
        flush=True,
    )
    print("=" * 78, flush=True)

    models, _, _, load_seconds = load_sample_models(sample)
    print(
        f"Models loaded: {len(models)} MAGs in {load_seconds:.3f}s.",
        flush=True,
    )

    model_sources = []
    for entry in sample["mags"]:
        model_path = Path(entry["model_path"]).expanduser().resolve()
        stat = model_path.stat()
        model_sources.append(
            {
                "mag_id": str(entry["mag_id"]),
                "path": str(model_path),
                "size_bytes": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        )

    medium_path = args.medium_file.expanduser().resolve()
    medium_stat = medium_path.stat()
    medium_fingerprint = {
        "path": str(medium_path),
        "size_bytes": int(medium_stat.st_size),
        "mtime_ns": int(medium_stat.st_mtime_ns),
    }

    build_start = perf_counter()
    community = LayeredCommunity(
        f"sample_{args.sample_index:03d}_product_scan",
        models,
        copy_models=False,
        create_solver=False,
    )
    build_seconds = perf_counter() - build_start
    print(
        f"Community built: {len(community.merged_model.reactions)} reactions; "
        f"{len(community.merged_model.metabolites)} metabolites; "
        f"{len(community.merged_model.genes)} genes in {build_seconds:.3f}s.",
        flush=True,
    )

    medium_start = perf_counter()
    medium, medium_audit, medium_file_audit = resolve_medium(
        community,
        medium_path,
        args.uptake_bound,
    )
    print(
        f"Medium resolved: {len(medium)} mapped exchanges in "
        f"{perf_counter() - medium_start:.3f}s.",
        flush=True,
    )

    all_exchange_reactions = sorted(get_community_exchanges(community))
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

    if not scan_candidates:
        raise RuntimeError("Product prefilter removed every global exchange.")

    scan_configuration = {
        "checkpoint_version": 2,
        "sample_id": str(sample["sample_id"]),
        "manifest_file": str(args.manifest.expanduser().resolve()),
        "model_sources": model_sources,
        "minimal_growth": float(args.minimal_growth),
        "production_tolerance": float(args.production_tolerance),
        "keep_oxygen": bool(args.keep_oxygen),
        "prefilter_before_optimization": True,
        "medium_file": str(medium_path),
        "medium_fingerprint": medium_fingerprint,
    }

    def scan_state_matches(payload: dict) -> bool:
        configuration = payload.get("configuration", {})
        return (
            int(configuration.get("checkpoint_version", -1)) == 2
            and float(configuration.get("minimal_growth", -1))
            == float(args.minimal_growth)
            and float(configuration.get("production_tolerance", -1))
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
            and str(configuration.get("sample_id", ""))
            == str(sample["sample_id"])
            and str(configuration.get("manifest_file", ""))
            == str(args.manifest.expanduser().resolve())
            and list(configuration.get("model_sources", []))
            == list(model_sources)
            and str(configuration.get("medium_file", ""))
            == str(medium_path)
            and dict(configuration.get("medium_fingerprint", {}))
            == dict(medium_fingerprint)
            and list(payload.get("candidate_exchange_reactions", []))
            == list(scan_candidates)
        )

    checkpoint_path = stage_dir / "product_scan.yaml"
    partial_path = stage_dir / "product_scan_partial.yaml"

    if checkpoint_path.exists():
        cached = yaml_load(checkpoint_path)
        if cached.get("status") == "complete" and scan_state_matches(cached):
            print(
                f"Product Scan checkpoint already complete: "
                f"{len(cached.get('exchange_results', {}))}/"
                f"{len(scan_candidates)} exchanges.",
                flush=True,
            )
            return

    completed_results = {}
    if partial_path.exists():
        partial = yaml_load(partial_path)
        if scan_state_matches(partial):
            completed_results = dict(partial.get("exchange_results", {}))

    lp_start = perf_counter()
    lp = build_product_scan_lp(
        community,
        minimal_growth=args.minimal_growth,
        medium=medium,
    )
    lp_build_seconds = perf_counter() - lp_start
    print(
        f"LP built: {len(lp.reaction_ids)} variables; "
        f"{lp.a_eq.shape[0]} equality rows; nnz={lp.a_eq.nnz}; "
        f"{len(scan_candidates)} product candidates in "
        f"{lp_build_seconds:.3f}s.",
        flush=True,
    )

    feasibility_start = perf_counter()
    feasibility = check_product_scan_feasibility(
        lp,
        tolerance=args.production_tolerance,
    )
    feasibility_seconds = perf_counter() - feasibility_start
    print(
        f"Feasibility check finished: optimal={feasibility['optimal']} "
        f"status={feasibility['status']} in {feasibility_seconds:.3f}s.",
        flush=True,
    )
    if not feasibility["optimal"]:
        raise RuntimeError(
            "Product Scan reference LP is infeasible/non-optimal: "
            + feasibility["message"]
        )

    exchange_results = dict(completed_results)

    def persist(rid: str, record: dict) -> None:
        exchange_results[rid] = record
        optimal = sum(
            1
            for value in exchange_results.values()
            if value.get("status") == "optimal"
            and value.get("maximum") is not None
        )
        yaml_dump_atomic(
            {
                "status": "partial",
                "configuration": scan_configuration,
                "candidate_exchange_reactions": list(scan_candidates),
                "number_total_exchange_reactions": len(all_exchange_reactions),
                "number_prefilter_candidates": len(scan_candidates),
                "number_prefilter_excluded": len(prefilter_excluded),
                "number_attempted": len(exchange_results),
                "number_completed_optimal": optimal,
                "number_remaining_optimal": len(scan_candidates) - optimal,
                "exchange_results": exchange_results,
            },
            partial_path,
        )

    unfiltered, exchange_results, scan_wall_seconds = run_parallel_product_scan(
        lp,
        scan_candidates,
        workers=args.workers,
        tolerance=args.production_tolerance,
        completed_results=exchange_results,
        progress_callback=persist,
    )

    nonoptimal = {
        rid: record
        for rid, record in exchange_results.items()
        if rid in scan_candidates
        and (
            record.get("status") != "optimal"
            or record.get("maximum") is None
        )
    }
    if nonoptimal:
        raise RuntimeError(
            "Product Scan contains non-optimal exchange solves: "
            + ", ".join(sorted(nonoptimal)[:20])
        )

    filtered, filter_audit = filter_product_candidates(
        community=community,
        max_secretion=unfiltered,
        keep_oxygen=args.keep_oxygen,
    )

    scan_seconds = sum(
        float(record.get("solve_seconds", 0.0))
        for rid, record in exchange_results.items()
        if rid in scan_candidates
    )

    payload = {
        "status": "complete",
        "configuration": scan_configuration,
        "solver": {
            "name": "scipy.optimize.linprog",
            "backend": "HiGHS",
            "parallel_workers": int(args.workers),
            "gurobi_used": False,
        },
        "candidate_exchange_reactions": list(scan_candidates),
        "number_total_exchange_reactions": len(all_exchange_reactions),
        "number_prefilter_candidates": len(scan_candidates),
        "number_prefilter_excluded": len(prefilter_excluded),
        "number_attempted": len(exchange_results),
        "number_completed_optimal": len(scan_candidates),
        "prefilter_audit": prefilter_audit,
        "exchange_results": exchange_results,
        "unfiltered_max_secretion": {
            rid: float(value) for rid, value in unfiltered.items()
        },
        "filtered_max_secretion": {
            rid: float(value) for rid, value in filtered.items()
        },
        "filter_audit": filter_audit,
        "scan_seconds": float(scan_seconds),
        "current_run_scan_wall_seconds": float(scan_wall_seconds),
    }

    yaml_dump_atomic(payload, checkpoint_path)
    yaml_dump_atomic(payload, partial_path)
    yaml_dump_atomic(
        {
            "status": "complete",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": str(sample["sample_id"]),
                "number_modeled_mags": int(sample["number_modeled_mags"]),
            },
            "solver": payload["solver"],
            "medium": {**medium_file_audit, **medium_audit},
            "number_candidates": len(scan_candidates),
            "number_producible": len(unfiltered),
            "number_filtered_products": len(filtered),
            "timing_seconds": {
                "model_loading": float(load_seconds),
                "community_build": float(build_seconds),
                "lp_build": float(lp_build_seconds),
                "feasibility": float(feasibility_seconds),
                "product_scan_wall": float(scan_wall_seconds),
                "sum_exchange_solve_seconds": float(scan_seconds),
                "total": float(perf_counter() - total_start),
            },
        },
        stage_dir / "product_scan_complete.yaml",
    )

    print(
        f"PRODUCT SCAN COMPLETE: {len(scan_candidates)} optimized; "
        f"{len(unfiltered)} producible; {len(filtered)} filtered; "
        f"wall={scan_wall_seconds:.3f}s; workers={args.workers}.",
        flush=True,
    )


if __name__ == "__main__":
    main()
