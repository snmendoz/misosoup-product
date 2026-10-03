"""Aggregate per-sample ocean experiment outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from common import yaml_dump_atomic, yaml_load


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def blank_summary(index: int, sample_id: str, n_mags: int, status: str) -> dict:
    return {
        "sample_index": index,
        "sample_id": sample_id,
        "status": status,
        "number_modeled_mags": n_mags,
        "number_unfiltered_producible_exchanges": "",
        "number_filtered_products": "",
        "k_star": "",
        "q_star": "",
        "number_selected_products": "",
        "selected_products": "",
        "pfba_objective": "",
        "stage_c_q_sum": "",
        "minimum_size": "",
        "minimum_fraction_of_original": "",
        "number_minimum_communities": "",
        "enumeration_complete": "",
        "minimum_communities": "",
        "individual_product_scan_seconds": "",
        "stage_a_seconds": "",
        "stage_b1_seconds": "",
        "stage_c_seconds": "",
        "stage_d_seconds": "",
        "total_seconds": "",
        "error_type": "",
        "error": "",
    }


def main() -> None:
    args = parse_args()
    manifest = yaml_load(args.manifest.expanduser().resolve())
    output_root = args.output_root.expanduser().resolve()
    samples = manifest["samples"]

    summary_rows = []
    producible_rows = []
    selected_rows = []
    community_rows = []
    member_rows = []
    failures = []

    for sample in samples:
        index = int(sample["index"])
        sample_id = str(sample["sample_id"])
        n_mags = int(sample["number_modeled_mags"])
        sample_dir = output_root / "samples" / f"{index:03d}"
        result_path = sample_dir / "result.yaml"
        failure_path = sample_dir / "failure.yaml"

        if result_path.exists():
            result = yaml_load(result_path)
            a = result["stage_a"]
            b = result["stage_b"]
            c = result["stage_c"]
            d = result["stage_d"]
            timing = result["timing_seconds"]
            nmin = int(d["minimum_size"])

            row = blank_summary(index, sample_id, n_mags, "ok")
            row.update(
                {
                    "number_unfiltered_producible_exchanges": result[
                        "individual_product_maxima"
                    ]["number_unfiltered_producible_exchanges"],
                    "number_filtered_products": result[
                        "individual_product_maxima"
                    ]["number_filtered_products"],
                    "k_star": int(a["k_star"]),
                    "q_star": float(b["q_star"]),
                    "number_selected_products": len(c["selected_products"]),
                    "selected_products": ";".join(c["selected_products"]),
                    "pfba_objective": float(c["pfba_objective"]),
                    "stage_c_q_sum": float(c["stage_c_q_sum"]),
                    "minimum_size": nmin,
                    "minimum_fraction_of_original": (
                        nmin / n_mags if n_mags else ""
                    ),
                    "number_minimum_communities": int(d["number_communities"]),
                    "enumeration_complete": bool(d["enumeration_complete"]),
                    "minimum_communities": ";".join(
                        "+".join(comm["organisms"])
                        for comm in d["communities"]
                    ),
                    "individual_product_scan_seconds": float(
                        timing["individual_product_scan"]
                    ),
                    "stage_a_seconds": float(timing["stage_a"]),
                    "stage_b1_seconds": float(timing["stage_b1"]),
                    "stage_c_seconds": float(timing["stage_c"]),
                    "stage_d_seconds": float(timing["stage_d"]),
                    "total_seconds": float(timing["total"]),
                }
            )
            summary_rows.append(row)

            for rid, maximum in result["individual_product_maxima"][
                "filtered_max_secretion"
            ].items():
                producible_rows.append(
                    {
                        "sample_index": index,
                        "sample_id": sample_id,
                        "product_id": rid,
                        "maximum_secretion": float(maximum),
                    }
                )

            for rid, info in c["reference_products"].items():
                selected_rows.append(
                    {
                        "sample_index": index,
                        "sample_id": sample_id,
                        "product_id": rid,
                        "max_individual": float(info["max_individual"]),
                        "required_fraction": float(info["required_fraction"]),
                        "required_flux": float(info["required_flux"]),
                        "reference_flux": float(info["reference_flux"]),
                        "normalized_flux": float(info["normalized_flux"]),
                        "q_value": float(info["q_value"]),
                    }
                )

            for ci, comm in enumerate(d["communities"], start=1):
                cid = f"{index:03d}_C{ci:04d}"
                community_rows.append(
                    {
                        "sample_index": index,
                        "sample_id": sample_id,
                        "community_index": ci,
                        "community_id": cid,
                        "size": int(comm["size"]),
                        "mags": ";".join(comm["organisms"]),
                        "community_growth": float(comm["community_growth"]),
                    }
                )
                for mag_id in comm["organisms"]:
                    member_rows.append(
                        {
                            "sample_index": index,
                            "sample_id": sample_id,
                            "community_index": ci,
                            "community_id": cid,
                            "mag_id": mag_id,
                            "organism_growth": float(
                                comm["organism_growth"][mag_id]
                            ),
                        }
                    )

        elif failure_path.exists():
            failure = yaml_load(failure_path)
            failures.append(failure)
            row = blank_summary(index, sample_id, n_mags, "failed")
            row["total_seconds"] = float(
                failure.get("timing_seconds", {}).get("total", 0.0)
            )
            row["error_type"] = failure.get("error_type", "")
            row["error"] = failure.get("error", "")
            summary_rows.append(row)

        else:
            missing = {
                "status": "missing",
                "sample": {"index": index, "sample_id": sample_id},
                "error_type": "MissingResult",
                "error": (
                    "No result.yaml or failure.yaml was found; "
                    "the array task may have been cancelled or timed out."
                ),
            }
            failures.append(missing)
            row = blank_summary(index, sample_id, n_mags, "missing")
            row["error_type"] = "MissingResult"
            row["error"] = missing["error"]
            summary_rows.append(row)

    summary_rows.sort(key=lambda x: int(x["sample_index"]))

    summary_fields = [
        "sample_index",
        "sample_id",
        "status",
        "number_modeled_mags",
        "number_unfiltered_producible_exchanges",
        "number_filtered_products",
        "k_star",
        "q_star",
        "number_selected_products",
        "selected_products",
        "pfba_objective",
        "stage_c_q_sum",
        "minimum_size",
        "minimum_fraction_of_original",
        "number_minimum_communities",
        "enumeration_complete",
        "minimum_communities",
        "individual_product_scan_seconds",
        "stage_a_seconds",
        "stage_b1_seconds",
        "stage_c_seconds",
        "stage_d_seconds",
        "total_seconds",
        "error_type",
        "error",
    ]

    write_csv(output_root / "sample_summary.csv", summary_fields, summary_rows)
    write_csv(
        output_root / "producible_products_long.csv",
        ["sample_index", "sample_id", "product_id", "maximum_secretion"],
        producible_rows,
    )
    write_csv(
        output_root / "selected_products_long.csv",
        [
            "sample_index",
            "sample_id",
            "product_id",
            "max_individual",
            "required_fraction",
            "required_flux",
            "reference_flux",
            "normalized_flux",
            "q_value",
        ],
        selected_rows,
    )
    write_csv(
        output_root / "minimal_communities.csv",
        [
            "sample_index",
            "sample_id",
            "community_index",
            "community_id",
            "size",
            "mags",
            "community_growth",
        ],
        community_rows,
    )
    write_csv(
        output_root / "minimal_community_members_long.csv",
        [
            "sample_index",
            "sample_id",
            "community_index",
            "community_id",
            "mag_id",
            "organism_growth",
        ],
        member_rows,
    )

    successful = sum(row["status"] == "ok" for row in summary_rows)
    failed = len(summary_rows) - successful
    aggregate = {
        "number_samples_expected": len(samples),
        "number_samples_successful": successful,
        "number_samples_failed_or_missing": failed,
        "failures": failures,
        "files": {
            "sample_summary": str(output_root / "sample_summary.csv"),
            "producible_products_long": str(
                output_root / "producible_products_long.csv"
            ),
            "selected_products_long": str(
                output_root / "selected_products_long.csv"
            ),
            "minimal_communities": str(
                output_root / "minimal_communities.csv"
            ),
            "minimal_community_members_long": str(
                output_root / "minimal_community_members_long.csv"
            ),
        },
    }
    yaml_dump_atomic(aggregate, output_root / "aggregate.yaml")

    print("=" * 72)
    print("Ocean-sample aggregation complete")
    print("=" * 72)
    print(f"Expected       : {len(samples)}")
    print(f"Successful     : {successful}")
    print(f"Failed/missing : {failed}")
    print(f"Summary        : {output_root / 'sample_summary.csv'}")
    print(f"Products       : {output_root / 'selected_products_long.csv'}")
    print(f"Min communities: {output_root / 'minimal_communities.csv'}")
    print("=" * 72)


if __name__ == "__main__":
    main()
