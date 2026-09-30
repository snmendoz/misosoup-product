"""Run individual-growth QC for a chunk of ocean MAG models."""

from __future__ import annotations

import argparse
import json
import math
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
from misosoup.library.readwrite import load_models
from misosoup.reframed.layered_community import LayeredCommunity


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--chunk-index", required=True, type=int)
    parser.add_argument("--chunk-size", required=True, type=int)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--medium-file", required=True, type=Path)
    parser.add_argument("--minimal-growth", type=float, default=0.01)
    parser.add_argument("--uptake-bound", type=float, default=-1000.0)
    parser.add_argument("--rich-uptake-bound", type=float, default=-1000.0)
    parser.add_argument("--production-tolerance", type=float, default=1e-6)
    parser.add_argument("--integer-tolerance", type=float, default=1e-9)
    return parser.parse_args()


def solver_params(prod_tol: float, int_tol: float) -> dict:
    return {
        Parameter.OPTIMALITY_TOL: prod_tol,
        Parameter.FEASIBILITY_TOL: prod_tol,
        Parameter.INT_FEASIBILITY_TOL: int_tol,
    }


def solve_max_growth(
    model,
    community_id: str,
    medium_entries: list[dict],
    params: dict,
    fallback_uptake_bound: float,
) -> dict:
    community = LayeredCommunity(
        community_id,
        [model],
        copy_models=False,
        params=params,
    )

    medium, audit = resolve_medium_for_community(
        community,
        medium_entries,
        fallback_uptake_bound,
    )

    community.setup_medium(medium)

    solution = community.solver.solve(
        objective={"community_growth": 1},
        get_values=["community_growth"],
        minimize=False,
    )

    growth = None

    if solution.status == Status.OPTIMAL:
        growth = float(
            solution.values.get(
                "community_growth",
                0.0,
            )
        )

    return {
        "status": str(solution.status),
        "max_growth": growth,
        "number_medium_tokens": int(
            audit["number_tokens"]
        ),
        "number_medium_tokens_matched": int(
            audit["number_matched_tokens"]
        ),
        "number_medium_tokens_unresolved": int(
            audit["number_unresolved_tokens"]
        ),
        "matched": audit["matched"],
        "unresolved": audit["unresolved"],
    }


def solve_rich_growth(
    model,
    community_id: str,
    params: dict,
    rich_uptake_bound: float,
) -> dict:
    community = LayeredCommunity(
        community_id,
        [model],
        copy_models=False,
        params=params,
    )

    global_exchanges = [
        rid
        for rid in community.merged_model.reactions
        if rid.startswith("R_EX_")
        and not rid.endswith("_i")
    ]

    rich_medium = {
        rid: float(rich_uptake_bound)
        for rid in global_exchanges
    }

    community.setup_medium(rich_medium)

    solution = community.solver.solve(
        objective={"community_growth": 1},
        get_values=["community_growth"],
        minimize=False,
    )

    growth = None

    if solution.status == Status.OPTIMAL:
        growth = float(
            solution.values.get(
                "community_growth",
                0.0,
            )
        )

    return {
        "status": str(solution.status),
        "max_growth": growth,
        "number_global_exchanges_opened": len(
            global_exchanges
        ),
        "uptake_bound": float(
            rich_uptake_bound
        ),
    }


def classify(
    selected: dict,
    rich: dict,
    minimal_growth: float,
    tolerance: float,
) -> str:
    selected_growth = selected.get("max_growth")
    rich_growth = rich.get("max_growth")

    selected_pass = (
        selected.get("status") == str(Status.OPTIMAL)
        and selected_growth is not None
        and selected_growth + tolerance >= minimal_growth
    )

    rich_pass = (
        rich.get("status") == str(Status.OPTIMAL)
        and rich_growth is not None
        and rich_growth + tolerance >= minimal_growth
    )

    if selected_pass:
        return "grows_in_medium"

    if rich_pass:
        return "medium_limited"

    if rich.get("status") == str(Status.OPTIMAL):
        return "no_growth_even_rich"

    return "solver_error"


def write_json_atomic(data: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")

    with tmp.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)

    tmp.replace(path)


def run_one(
    entry: dict,
    output_root: Path,
    medium_entries: list[dict],
    medium_audit: dict,
    params: dict,
    minimal_growth: float,
    uptake_bound: float,
    rich_uptake_bound: float,
    tolerance: float,
) -> dict:
    index = int(entry["index"])
    mag_id = str(entry["mag_id"])
    model_path = Path(entry["model_path"]).expanduser().resolve()
    model_dir = output_root / "models" / f"{index:04d}"
    model_dir.mkdir(parents=True, exist_ok=True)

    start = perf_counter()

    try:
        models = load_models([str(model_path)])

        if len(models) != 1:
            raise RuntimeError(
                f"Expected one model, loaded {len(models)}"
            )

        model = models[0]
        original_model_id = str(model.id)
        biomass_reaction = str(model.biomass_reaction)
        model.id = "qc_org"

        selected = solve_max_growth(
            model=model,
            community_id=f"qc_selected_{index:04d}",
            medium_entries=medium_entries,
            params=params,
            fallback_uptake_bound=uptake_bound,
        )

        rich = solve_rich_growth(
            model=model,
            community_id=f"qc_rich_{index:04d}",
            params=params,
            rich_uptake_bound=rich_uptake_bound,
        )

        selected_growth = selected["max_growth"]
        rich_growth = rich["max_growth"]

        selected_pass = (
            selected["status"] == str(Status.OPTIMAL)
            and selected_growth is not None
            and selected_growth + tolerance >= minimal_growth
        )

        rich_pass = (
            rich["status"] == str(Status.OPTIMAL)
            and rich_growth is not None
            and rich_growth + tolerance >= minimal_growth
        )

        diagnostic_class = classify(
            selected=selected,
            rich=rich,
            minimal_growth=minimal_growth,
            tolerance=tolerance,
        )

        result = {
            "status": "ok",
            "index": index,
            "mag_id": mag_id,
            "model_path": str(model_path),
            "original_model_id": original_model_id,
            "biomass_reaction": biomass_reaction,
            "minimal_growth_requirement": float(
                minimal_growth
            ),
            "medium_file": medium_audit,
            "selected_medium": selected,
            "rich_medium": rich,
            "grows_in_selected_medium": bool(
                selected_pass
            ),
            "grows_in_rich_medium": bool(
                rich_pass
            ),
            "diagnostic_class": diagnostic_class,
            "elapsed_seconds": perf_counter() - start,
        }

        yaml_dump_atomic(
            result,
            model_dir / "result.yaml",
        )

        write_json_atomic(
            {
                "status": "ok",
                "mag_id": mag_id,
                "diagnostic_class": diagnostic_class,
                "selected_medium_max_growth": selected_growth,
                "rich_medium_max_growth": rich_growth,
            },
            model_dir / "status.json",
        )

        print(
            f"[{index}] {mag_id}: "
            f"selected={selected_growth} "
            f"rich={rich_growth} "
            f"class={diagnostic_class}",
            flush=True,
        )

        return result

    except Exception as error:
        trace = traceback.format_exc()

        failure = {
            "status": "failed",
            "index": index,
            "mag_id": mag_id,
            "model_path": str(model_path),
            "error_type": type(error).__name__,
            "error": str(error),
            "traceback": trace,
            "elapsed_seconds": perf_counter() - start,
        }

        yaml_dump_atomic(
            failure,
            model_dir / "failure.yaml",
        )

        write_json_atomic(
            {
                "status": "failed",
                "mag_id": mag_id,
                "error_type": type(error).__name__,
                "error": str(error),
            },
            model_dir / "status.json",
        )

        print(
            f"[{index}] {mag_id}: FAILED "
            f"{type(error).__name__}: {error}",
            flush=True,
        )

        return failure


def main() -> None:
    args = parse_args()

    if args.chunk_size < 1:
        raise ValueError("--chunk-size must be >= 1")

    if args.minimal_growth <= 0:
        raise ValueError("--minimal-growth must be > 0")

    if args.uptake_bound >= 0:
        raise ValueError("--uptake-bound must be negative")

    if args.rich_uptake_bound >= 0:
        raise ValueError("--rich-uptake-bound must be negative")

    manifest = yaml_load(
        args.manifest.expanduser().resolve()
    )

    models = manifest["models"]
    start_index = args.chunk_index * args.chunk_size
    stop_index = min(
        start_index + args.chunk_size,
        len(models),
    )

    if start_index >= len(models):
        raise IndexError(
            f"Chunk {args.chunk_index} starts at {start_index}, "
            f"but manifest has only {len(models)} models."
        )

    output_root = args.output_root.expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    medium_entries, medium_audit = read_medium_spec(
        args.medium_file,
        default_uptake_bound=args.uptake_bound,
    )

    params = solver_params(
        args.production_tolerance,
        args.integer_tolerance,
    )

    print("=" * 72)
    print(
        f"Model QC chunk {args.chunk_index}: "
        f"models {start_index}..{stop_index - 1}"
    )
    print(
        f"Medium: {args.medium_file.expanduser().resolve()}"
    )
    print("=" * 72, flush=True)

    results = []

    for entry in models[start_index:stop_index]:
        results.append(
            run_one(
                entry=entry,
                output_root=output_root,
                medium_entries=medium_entries,
                medium_audit=medium_audit,
                params=params,
                minimal_growth=args.minimal_growth,
                uptake_bound=args.uptake_bound,
                rich_uptake_bound=args.rich_uptake_bound,
                tolerance=args.production_tolerance,
            )
        )

    chunk_dir = output_root / "chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    yaml_dump_atomic(
        {
            "chunk_index": int(args.chunk_index),
            "start_index": int(start_index),
            "stop_index_exclusive": int(stop_index),
            "number_models": int(len(results)),
            "number_ok": int(
                sum(
                    result["status"] == "ok"
                    for result in results
                )
            ),
            "number_failed": int(
                sum(
                    result["status"] != "ok"
                    for result in results
                )
            ),
        },
        chunk_dir / f"chunk_{args.chunk_index:04d}.yaml",
    )


if __name__ == "__main__":
    main()
