"""Stage 3: minimum product-preserving communities plus final audits."""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from common import yaml_dump_atomic, yaml_load
from misosoup.library.minimal_product_communities import findMinimalProductCommunities
from misosoup.reframed.layered_community import LayeredCommunity
from run_sample import (
    build_medium_uptake_audit,
    build_product_preservation_audit,
    map_stage_d,
    write_csv_atomic,
)
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
    parser.add_argument("--product-retention", type=float, default=0.90)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument("--max-minimal-communities", type=int, default=100)
    return parser.parse_args()


def main():
    args = parse_args()
    _, sample = get_sample(args.manifest, args.sample_index)
    root = sample_dir(args.output_root, args.sample_index)
    stage_dir = root / "03_stage_d"
    stage_dir.mkdir(parents=True, exist_ok=True)
    total_start = perf_counter()

    try:
        ab_path = root / "01_reference_ab" / "reference_ab.yaml"
        c_path = root / "02_stage_c" / "stage_c.yaml"

        if not ab_path.exists():
            raise FileNotFoundError(f"Missing Stage A/B checkpoint: {ab_path}")
        if not c_path.exists():
            raise FileNotFoundError(f"Missing Stage C checkpoint: {c_path}")

        ab = yaml_load(ab_path)
        c = yaml_load(c_path)

        if ab.get("status") != "complete" or c.get("status") != "complete":
            raise RuntimeError("Upstream staged checkpoint is not complete.")

        models, solver_to_mag, _, load_seconds = load_sample_models(sample)
        params = solver_params(args.production_tolerance, args.integer_tolerance)

        t = perf_counter()
        community = LayeredCommunity(
            f"sample_{args.sample_index:03d}_stage_d",
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
        stage_d_raw = findMinimalProductCommunities(
            community=community,
            medium=medium,
            reference_products=c["stage_c"]["reference_products"],
            minimal_growth=args.minimal_growth,
            product_retention=args.product_retention,
            tolerance=args.production_tolerance,
            max_communities=args.max_minimal_communities,
        )
        stage_d_seconds = perf_counter() - t
        stage_d = map_stage_d(stage_d_raw, solver_to_mag)

        product_rows = build_product_preservation_audit(
            stage_d=stage_d,
            stage_c=c["stage_c"],
            full_product_maxima=ab[
                "individual_product_maxima"
            ]["filtered_max_secretion"],
            community=community,
            tolerance=args.production_tolerance,
        )
        medium_rows = build_medium_uptake_audit(
            stage_d=stage_d,
            medium=medium,
            medium_audit=medium_audit,
            community=community,
            tolerance=args.production_tolerance,
        )

        write_csv_atomic(
            product_rows,
            [
                "community_index",
                "community_size",
                "members",
                "product_reaction",
                "metabolite_id",
                "metabolite_name",
                "metabolite_formula",
                "max_flux_full_community",
                "stage_c_reference_flux",
                "required_retention_fraction",
                "stage_d_required_flux",
                "stage_d_observed_flux",
                "observed_retention_fraction",
                "margin_to_requirement",
                "passes_requirement",
            ],
            root / "product_preservation_audit.csv",
        )

        write_csv_atomic(
            medium_rows,
            [
                "community_index",
                "community_size",
                "members",
                "medium_tokens",
                "medium_names",
                "exchange_reaction",
                "metabolite_id",
                "metabolite_name",
                "metabolite_formula",
                "allowed_lower_bound",
                "allowed_uptake_magnitude",
                "exchange_flux",
                "actual_uptake_magnitude",
                "actual_secretion_flux",
                "fraction_of_allowed_uptake",
                "margin_from_lower_bound",
                "at_uptake_limit",
            ],
            root / "medium_uptake_audit.csv",
        )

        ab_t = ab["timing_seconds"]
        c_t = c["timing_seconds"]
        stage_d_total = perf_counter() - total_start
        total_compute = float(ab_t["total"]) + float(c_t["total"]) + stage_d_total

        result = {
            "status": "ok",
            "sample": {
                "index": int(args.sample_index),
                "sample_id": str(sample["sample_id"]),
                "number_present_mags_in_matrix": int(
                    sample["number_present_mags_in_matrix"]
                ),
                "number_modeled_mags": int(sample["number_modeled_mags"]),
                "members": sample["mags"],
            },
            "configuration": {
                **ab["configuration"],
                "product_retention": float(args.product_retention),
                "max_minimal_communities": int(args.max_minimal_communities),
            },
            "medium": {
                **medium_file_audit,
                **medium_audit,
            },
            "individual_product_maxima": ab["individual_product_maxima"],
            "stage_a": ab["stage_a"],
            "stage_b": ab["stage_b"],
            "stage_c": c["stage_c"],
            "stage_d": stage_d,
            "audits": {
                "product_preservation_rows": len(product_rows),
                "medium_uptake_rows": len(medium_rows),
                "product_preservation_file": str(
                    root / "product_preservation_audit.csv"
                ),
                "medium_uptake_file": str(root / "medium_uptake_audit.csv"),
            },
            "timing_seconds": {
                "individual_product_scan": float(
                    ab_t["individual_product_scan"]
                ),
                "stage_a": float(ab_t["stage_a"]),
                "stage_b1": float(ab_t["stage_b1"]),
                "stage_c": float(c_t["stage_c"]),
                "stage_d": float(stage_d_seconds),
                "stage_ab_total": float(ab_t["total"]),
                "stage_c_total": float(c_t["total"]),
                "stage_d_total": float(stage_d_total),
                "total": float(total_compute),
            },
        }

        yaml_dump_atomic(result, root / "result.yaml")
        yaml_dump_atomic(
            {
                "status": "complete",
                "stage_d": stage_d,
                "timing_seconds": {
                    "model_loading": float(load_seconds),
                    "community_build": float(build_seconds),
                    "stage_d": float(stage_d_seconds),
                    "total": float(stage_d_total),
                },
            },
            stage_dir / "stage_d.yaml",
        )
        mark_complete(stage_dir)

        failure_path = root / "failure.yaml"
        if failure_path.exists():
            failure_path.unlink()

        print(
            f"STAGE D COMPLETE in {stage_d_total:.3f}s; "
            f"N_min={stage_d['minimum_size']}; "
            f"communities={stage_d['number_communities']}.",
            flush=True,
        )

    except Exception as error:
        write_stage_failure(
            stage_name="stage_d",
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
