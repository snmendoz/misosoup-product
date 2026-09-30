"""Collect individual ocean-model QC results into CSV and YAML summaries."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from common import yaml_dump_atomic, yaml_load


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    manifest = yaml_load(
        args.manifest.expanduser().resolve()
    )

    output_root = args.output_root.expanduser().resolve()
    rows = []
    classes = Counter()

    for entry in manifest["models"]:
        index = int(entry["index"])
        mag_id = str(entry["mag_id"])
        model_dir = output_root / "models" / f"{index:04d}"
        result_path = model_dir / "result.yaml"
        failure_path = model_dir / "failure.yaml"

        if result_path.exists():
            result = yaml_load(result_path)
            selected = result["selected_medium"]
            rich = result["rich_medium"]
            diagnostic_class = str(
                result["diagnostic_class"]
            )
            classes[diagnostic_class] += 1

            rows.append(
                {
                    "index": index,
                    "mag_id": mag_id,
                    "status": "ok",
                    "diagnostic_class": diagnostic_class,
                    "biomass_reaction": result[
                        "biomass_reaction"
                    ],
                    "selected_medium_max_growth": selected[
                        "max_growth"
                    ],
                    "grows_in_selected_medium": result[
                        "grows_in_selected_medium"
                    ],
                    "medium_tokens_matched": selected[
                        "number_medium_tokens_matched"
                    ],
                    "medium_tokens_total": selected[
                        "number_medium_tokens"
                    ],
                    "rich_medium_max_growth": rich[
                        "max_growth"
                    ],
                    "grows_in_rich_medium": result[
                        "grows_in_rich_medium"
                    ],
                    "rich_global_exchanges_opened": rich[
                        "number_global_exchanges_opened"
                    ],
                    "elapsed_seconds": result[
                        "elapsed_seconds"
                    ],
                    "error_type": "",
                    "error": "",
                    "model_path": result["model_path"],
                }
            )

        elif failure_path.exists():
            failure = yaml_load(failure_path)
            classes["solver_error"] += 1

            rows.append(
                {
                    "index": index,
                    "mag_id": mag_id,
                    "status": "failed",
                    "diagnostic_class": "solver_error",
                    "biomass_reaction": "",
                    "selected_medium_max_growth": "",
                    "grows_in_selected_medium": "",
                    "medium_tokens_matched": "",
                    "medium_tokens_total": "",
                    "rich_medium_max_growth": "",
                    "grows_in_rich_medium": "",
                    "rich_global_exchanges_opened": "",
                    "elapsed_seconds": failure.get(
                        "elapsed_seconds",
                        "",
                    ),
                    "error_type": failure.get(
                        "error_type",
                        "",
                    ),
                    "error": failure.get(
                        "error",
                        "",
                    ),
                    "model_path": failure.get(
                        "model_path",
                        entry["model_path"],
                    ),
                }
            )

        else:
            classes["missing"] += 1

            rows.append(
                {
                    "index": index,
                    "mag_id": mag_id,
                    "status": "missing",
                    "diagnostic_class": "missing",
                    "biomass_reaction": "",
                    "selected_medium_max_growth": "",
                    "grows_in_selected_medium": "",
                    "medium_tokens_matched": "",
                    "medium_tokens_total": "",
                    "rich_medium_max_growth": "",
                    "grows_in_rich_medium": "",
                    "rich_global_exchanges_opened": "",
                    "elapsed_seconds": "",
                    "error_type": "",
                    "error": "",
                    "model_path": entry["model_path"],
                }
            )

    rows.sort(key=lambda row: int(row["index"]))

    csv_path = output_root / "model_qc.csv"

    fieldnames = [
        "index",
        "mag_id",
        "status",
        "diagnostic_class",
        "biomass_reaction",
        "selected_medium_max_growth",
        "grows_in_selected_medium",
        "medium_tokens_matched",
        "medium_tokens_total",
        "rich_medium_max_growth",
        "grows_in_rich_medium",
        "rich_global_exchanges_opened",
        "elapsed_seconds",
        "error_type",
        "error",
        "model_path",
    ]

    with csv_path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(rows)

    aggregate = {
        "number_models": len(rows),
        "counts_by_diagnostic_class": dict(
            sorted(classes.items())
        ),
        "number_ok": sum(
            row["status"] == "ok"
            for row in rows
        ),
        "number_failed": sum(
            row["status"] == "failed"
            for row in rows
        ),
        "number_missing": sum(
            row["status"] == "missing"
            for row in rows
        ),
        "model_qc_csv": str(csv_path),
    }

    yaml_dump_atomic(
        aggregate,
        output_root / "aggregate.yaml",
    )

    print("=" * 72)
    print("Individual model QC collected")
    print("=" * 72)
    print(f"Models        : {len(rows)}")
    print(f"CSV           : {csv_path}")
    print("Classes:")
    for key, value in sorted(classes.items()):
        print(f"  {key}: {value}")
    print("=" * 72, flush=True)


if __name__ == "__main__":
    main()
