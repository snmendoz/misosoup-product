"""Run low-memory growth diagnostics once for every unique MAG in a manifest."""

from __future__ import annotations

import argparse
from pathlib import Path
from time import perf_counter

from misosoup.cobra.solver import Status

from common import (
    read_medium_spec,
    resolve_medium_for_community,
    yaml_dump_atomic,
    yaml_load,
)
from misosoup.library.readwrite import load_models
from misosoup.cobra.layered_community import LayeredCommunity
from staged_common import solver_params


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    return parser.parse_args()


def main():
    args = parse_args()
    manifest = yaml_load(args.manifest.expanduser().resolve())
    output_root = args.output_root.expanduser().resolve()
    qc_dir = output_root / "00_global_qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    unique = {}
    for sample in manifest["samples"]:
        for entry in sample["mags"]:
            mag_id = str(entry["mag_id"])
            model_path = str(Path(entry["model_path"]).expanduser().resolve())
            previous = unique.get(mag_id)
            if previous is not None and previous != model_path:
                raise RuntimeError(
                    f"MAG {mag_id} maps to multiple model files: "
                    f"{previous} and {model_path}"
                )
            unique[mag_id] = model_path

    medium_entries, medium_file_audit = read_medium_spec(
        args.medium_file.expanduser().resolve(),
        default_uptake_bound=args.uptake_bound,
    )
    params = solver_params(args.production_tolerance, args.integer_tolerance)

    print(
        f"Global QC: {len(unique)} unique MAG models; "
        f"growth threshold={args.minimal_growth}.",
        flush=True,
    )

    start = perf_counter()
    rows = []

    for index, (mag_id, model_path) in enumerate(sorted(unique.items()), start=1):
        entry = {
            "mag_id": mag_id,
            "model_path": model_path,
            "minimal_growth_requirement": float(args.minimal_growth),
        }

        try:
            model = load_models([model_path])[0]
            model.id = "org_qc"
            community = LayeredCommunity(
                f"qc_{index:05d}",
                [model],
                copy_models=False,
                params=params,
            )
            medium, audit = resolve_medium_for_community(
                community,
                medium_entries,
                args.uptake_bound,
            )
            community.setup_medium(medium)
            solution = community.solver.solve(
                objective={"community_growth": 1},
                get_values=["community_growth"],
                minimize=False,
            )
            entry["status"] = str(solution.status)
            entry["number_medium_tokens_matched"] = int(
                audit["number_matched_tokens"]
            )
            if solution.status == Status.OPTIMAL:
                growth = float(solution.values.get("community_growth", 0.0))
                entry["max_growth"] = growth
                entry["passes_minimal_growth"] = bool(
                    growth + args.production_tolerance >= args.minimal_growth
                )
            else:
                entry["max_growth"] = None
                entry["passes_minimal_growth"] = False
        except Exception as error:
            entry["status"] = "diagnostic_failed"
            entry["max_growth"] = None
            entry["passes_minimal_growth"] = False
            entry["error_type"] = type(error).__name__
            entry["error"] = str(error)

        rows.append(entry)

        if index % 25 == 0 or index == len(unique):
            passed = sum(bool(row["passes_minimal_growth"]) for row in rows)
            print(
                f"Global QC progress: {index}/{len(unique)}; "
                f"passing={passed}; elapsed={perf_counter() - start:.1f}s.",
                flush=True,
            )
            yaml_dump_atomic(
                {
                    "status": "running" if index < len(unique) else "complete",
                    "number_models": len(unique),
                    "number_processed": index,
                    "number_passing": passed,
                    "number_not_passing": index - passed,
                    "medium": medium_file_audit,
                    "models": rows,
                    "elapsed_seconds": float(perf_counter() - start),
                },
                qc_dir / "global_qc.yaml",
            )

    (qc_dir / "COMPLETE").write_text("ok\n", encoding="utf-8")


if __name__ == "__main__":
    main()
