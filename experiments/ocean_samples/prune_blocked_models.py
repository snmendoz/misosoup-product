"""Remove structurally blocked reactions from each unique MAG model.

Blocked reactions are identified by FVA with every environmental exchange
opened to +/- open_bound. The biomass reaction is deliberately excluded from
exchange opening because gapseq can encode biomass with an R_EX_* identifier;
instead, biomass is explicitly constrained to [0, open_bound].

The reduced models are written as:
    <model_stem>_without_blocked_reactions.xml

A transformed manifest pointing to the reduced models is produced only after
all unique models finish successfully. Per-model audit YAML files make the
preprocessing resumable.
"""

from __future__ import annotations

import argparse
import copy
import csv
import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from time import perf_counter

from gurobipy import Env
from reframed import FBA, FVA, ReactionType, save_cbmodel
from reframed.solvers.solution import Status

from common import yaml_dump_atomic, yaml_load
from misosoup.library.readwrite import load_models
from misosoup.reframed.gurobi_env_solver import GurobiEnvSolver


_WORKER_ENV = None


def _worker_init() -> None:
    """Create exactly one persistent Gurobi WLS environment per worker."""
    global _WORKER_ENV
    _WORKER_ENV = Env(
        params={
            "LogToConsole": 0,
            "Method": 1,
        }
    )


def _source_fingerprint(path: Path) -> dict:
    stat = path.stat()
    return {
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def _environmental_exchange_ids(model) -> list[str]:
    """Return environmental exchanges, explicitly excluding biomass."""
    biomass = model.biomass_reaction
    ids = []

    for rid, reaction in model.reactions.items():
        if rid == biomass:
            continue

        if (
            reaction.reaction_type == ReactionType.EXCHANGE
            or rid.startswith("R_EX")
        ):
            ids.append(rid)

    return sorted(set(ids))


def _open_exchange_constraints(model, open_bound: float) -> dict:
    """Open environmental exchanges and set biomass explicitly to [0, bound]."""
    constraints = {
        rid: (-float(open_bound), float(open_bound))
        for rid in _environmental_exchange_ids(model)
    }
    constraints[model.biomass_reaction] = (
        0.0,
        float(open_bound),
    )
    return constraints


def _dispose_solver(solver) -> None:
    try:
        solver.problem.dispose()
    except Exception:
        pass


def _max_biomass(model, constraints: dict) -> float:
    solver = GurobiEnvSolver(
        model=model,
        env=_WORKER_ENV,
    )
    try:
        solution = FBA(
            model,
            objective=model.biomass_reaction,
            constraints=constraints,
            solver=solver,
        )
    finally:
        _dispose_solver(solver)

    if solution.status != Status.OPTIMAL:
        raise RuntimeError(
            "Open-exchange biomass FBA is not optimal: "
            f"{solution.status}"
        )

    return float(solution.fobj)


def _remove_orphan_metabolites(model) -> list[str]:
    used = {
        metabolite_id
        for reaction in model.reactions.values()
        for metabolite_id in reaction.stoichiometry
    }
    orphan_ids = sorted(
        set(model.metabolites.keys()) - used
    )

    if orphan_ids:
        model.remove_metabolites(orphan_ids)

    return orphan_ids


def _audit_path(output_model: Path) -> Path:
    return output_model.with_name(
        output_model.stem + "_blocked_reactions.yaml"
    )


def _cache_valid(
    source_path: Path,
    output_model: Path,
    audit_path: Path,
    *,
    tolerance: float,
    open_bound: float,
) -> bool:
    if not output_model.exists() or not audit_path.exists():
        return False

    try:
        audit = yaml_load(audit_path)
    except Exception:
        return False

    return (
        audit.get("status") == "complete"
        and audit.get("source_model") == str(source_path)
        and audit.get("source_fingerprint") == _source_fingerprint(source_path)
        and float(audit.get("blocked_tolerance", -1.0))
        == float(tolerance)
        and float(audit.get("open_exchange_bound", -1.0))
        == float(open_bound)
        and audit.get("output_model") == str(output_model)
    )


def _prune_one(task: dict) -> dict:
    source_path = Path(task["source_path"]).resolve()
    output_model = Path(task["output_model"]).resolve()
    audit_path = _audit_path(output_model)
    tolerance = float(task["tolerance"])
    open_bound = float(task["open_bound"])

    output_model.parent.mkdir(parents=True, exist_ok=True)

    if _cache_valid(
        source_path,
        output_model,
        audit_path,
        tolerance=tolerance,
        open_bound=open_bound,
    ):
        cached = yaml_load(audit_path)
        cached["cache_reused"] = True
        return cached

    start = perf_counter()
    model = load_models([str(source_path)])[0]

    biomass = model.biomass_reaction
    if biomass not in model.reactions:
        raise RuntimeError(
            f"Biomass reaction {biomass!r} is absent from {source_path}."
        )

    original_reactions = len(model.reactions)
    original_metabolites = len(model.metabolites)
    exchange_ids = _environmental_exchange_ids(model)
    constraints = _open_exchange_constraints(model, open_bound)

    # Validate that the maximally permissive single-organism model is viable.
    biomass_before = _max_biomass(model, constraints)

    solver = GurobiEnvSolver(
        model=model,
        env=_WORKER_ENV,
    )
    try:
        variability = FVA(
            model,
            obj_frac=0,
            reactions=list(model.reactions.keys()),
            constraints=constraints,
            loopless=False,
            solver=solver,
        )
    finally:
        _dispose_solver(solver)

    blocked = sorted(
        rid
        for rid, bounds in variability.items()
        if (
            bounds[0] is not None
            and bounds[1] is not None
            and abs(float(bounds[0])) + abs(float(bounds[1]))
            < tolerance
        )
    )

    if biomass in blocked:
        raise RuntimeError(
            f"Biomass reaction {biomass} was classified as blocked "
            "with all environmental exchanges open."
        )

    for rid in blocked:
        model.remove_reaction(rid)

    orphan_metabolites = _remove_orphan_metabolites(model)

    remaining_constraints = {
        rid: bounds
        for rid, bounds in constraints.items()
        if rid in model.reactions
    }
    biomass_after = _max_biomass(
        model,
        remaining_constraints,
    )

    comparison_tolerance = max(
        tolerance,
        1e-7 * max(1.0, abs(biomass_before)),
    )
    if abs(biomass_after - biomass_before) > comparison_tolerance:
        raise RuntimeError(
            "Biomass validation changed after blocked-reaction pruning: "
            f"before={biomass_before:.12g}, after={biomass_after:.12g}, "
            f"allowed_difference={comparison_tolerance:.12g}."
        )

    temporary = output_model.with_name(
        output_model.stem + ".tmp.xml"
    )
    save_cbmodel(
        model,
        str(temporary),
        flavor="fbc2",
    )
    os.replace(temporary, output_model)

    audit = {
        "status": "complete",
        "cache_reused": False,
        "source_model": str(source_path),
        "source_fingerprint": _source_fingerprint(source_path),
        "output_model": str(output_model),
        "model_id": str(model.id),
        "biomass_reaction": str(biomass),
        "blocked_tolerance": tolerance,
        "open_exchange_bound": open_bound,
        "biomass_lower_bound_fva": 0.0,
        "biomass_upper_bound_fva": open_bound,
        "number_environmental_exchanges_opened": len(exchange_ids),
        "original_reactions": int(original_reactions),
        "blocked_reactions": len(blocked),
        "remaining_reactions": int(len(model.reactions)),
        "reaction_reduction_fraction": (
            len(blocked) / original_reactions
            if original_reactions
            else 0.0
        ),
        "original_metabolites": int(original_metabolites),
        "orphan_metabolites_removed": len(orphan_metabolites),
        "remaining_metabolites": int(len(model.metabolites)),
        "biomass_max_open_exchanges_before": biomass_before,
        "biomass_max_open_exchanges_after": biomass_after,
        "blocked_reaction_ids": blocked,
        "orphan_metabolite_ids": orphan_metabolites,
        "elapsed_seconds": float(perf_counter() - start),
    }
    yaml_dump_atomic(audit, audit_path)
    return audit


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output-model-dir", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--summary-csv", type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--blocked-tolerance", type=float, default=1e-9)
    parser.add_argument("--open-exchange-bound", type=float, default=1000.0)
    parser.add_argument(
        "--max-models",
        type=int,
        default=None,
        help="Pilot/debug limit. Omit for every unique model.",
    )
    return parser.parse_args()


def _write_summary(path: Path, audits: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source_model",
        "output_model",
        "model_id",
        "biomass_reaction",
        "original_reactions",
        "blocked_reactions",
        "remaining_reactions",
        "reaction_reduction_fraction",
        "original_metabolites",
        "orphan_metabolites_removed",
        "remaining_metabolites",
        "number_environmental_exchanges_opened",
        "biomass_lower_bound_fva",
        "biomass_upper_bound_fva",
        "biomass_max_open_exchanges_before",
        "biomass_max_open_exchanges_after",
        "elapsed_seconds",
        "cache_reused",
    ]

    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for audit in sorted(
            audits,
            key=lambda row: row["source_model"],
        ):
            writer.writerow(
                {
                    key: audit.get(key, "")
                    for key in fields
                }
            )

    os.replace(temporary, path)


def main() -> None:
    args = parse_args()

    if args.workers < 1 or args.workers > 2:
        raise ValueError(
            "--workers must be 1 or 2 because the current Gurobi WLS "
            "baseline supports at most two concurrent sessions."
        )
    if args.blocked_tolerance <= 0:
        raise ValueError("--blocked-tolerance must be > 0.")
    if args.open_exchange_bound <= 0:
        raise ValueError("--open-exchange-bound must be > 0.")

    manifest_path = args.manifest.expanduser().resolve()
    output_model_dir = args.output_model_dir.expanduser().resolve()
    output_manifest = args.output_manifest.expanduser().resolve()
    summary_csv = (
        args.summary_csv.expanduser().resolve()
        if args.summary_csv is not None
        else output_model_dir / "blocked_reaction_pruning_summary.csv"
    )

    manifest = yaml_load(manifest_path)

    model_to_mags: dict[str, set[str]] = {}
    for sample in manifest["samples"]:
        for entry in sample["mags"]:
            source = str(
                Path(entry["model_path"]).expanduser().resolve()
            )
            model_to_mags.setdefault(source, set()).add(
                str(entry["mag_id"])
            )

    source_paths = sorted(model_to_mags)
    if args.max_models is not None:
        if args.max_models < 1:
            raise ValueError("--max-models must be >= 1.")
        source_paths = source_paths[: args.max_models]

    print(
        f"Blocked-reaction pruning: {len(source_paths)} unique models; "
        f"workers={args.workers}; all environmental exchanges opened to "
        f"+/-{args.open_exchange_bound:g}; tolerance="
        f"{args.blocked_tolerance:g}.",
        flush=True,
    )

    tasks = []
    mapping = {}
    claimed_outputs = {}

    for source_text in source_paths:
        source = Path(source_text)
        output_model = output_model_dir / (
            source.stem + "_without_blocked_reactions.xml"
        )
        output_key = str(output_model.resolve())

        previous_source = claimed_outputs.get(output_key)
        if previous_source is not None and previous_source != source_text:
            raise RuntimeError(
                "Two source models would produce the same reduced filename: "
                f"{previous_source} and {source_text} -> {output_key}"
            )

        claimed_outputs[output_key] = source_text
        mapping[source_text] = output_key
        tasks.append(
            {
                "source_path": source_text,
                "output_model": str(output_model),
                "tolerance": args.blocked_tolerance,
                "open_bound": args.open_exchange_bound,
            }
        )

    start = perf_counter()
    audits = []

    context = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=args.workers,
        mp_context=context,
        initializer=_worker_init,
    ) as executor:
        futures = {
            executor.submit(_prune_one, task): task
            for task in tasks
        }

        for completed, future in enumerate(
            as_completed(futures),
            start=1,
        ):
            task = futures[future]
            try:
                audit = future.result()
            except Exception as error:
                raise RuntimeError(
                    "Blocked-reaction pruning failed for "
                    f"{task['source_path']}: {error}"
                ) from error

            audits.append(audit)
            print(
                f"Pruning {completed}/{len(tasks)}: "
                f"{Path(audit['source_model']).name} "
                f"{audit['original_reactions']} -> "
                f"{audit['remaining_reactions']} reactions "
                f"(-{audit['blocked_reactions']}); "
                f"{audit['original_metabolites']} -> "
                f"{audit['remaining_metabolites']} metabolites; "
                f"{audit['elapsed_seconds']:.2f}s"
                + (
                    " [cached]"
                    if audit.get("cache_reused")
                    else ""
                ),
                flush=True,
            )

            _write_summary(summary_csv, audits)

    if args.max_models is not None:
        print(
            "Pilot mode: --max-models was supplied, so no transformed "
            "manifest will be written.",
            flush=True,
        )
        return

    transformed = copy.deepcopy(manifest)
    for sample in transformed["samples"]:
        for entry in sample["mags"]:
            source = str(
                Path(entry["model_path"]).expanduser().resolve()
            )
            entry["original_model_path"] = source
            entry["model_path"] = mapping[source]

    transformed["blocked_reaction_pruning"] = {
        "source_manifest": str(manifest_path),
        "output_model_dir": str(output_model_dir),
        "blocked_tolerance": float(args.blocked_tolerance),
        "open_exchange_bound": float(args.open_exchange_bound),
        "number_unique_models": len(source_paths),
        "summary_csv": str(summary_csv),
        "elapsed_seconds": float(perf_counter() - start),
    }

    yaml_dump_atomic(transformed, output_manifest)

    print(
        f"Blocked-reaction pruning COMPLETE: "
        f"{len(source_paths)} models in {perf_counter() - start:.1f}s.",
        flush=True,
    )
    print(
        f"Reduced manifest: {output_manifest}",
        flush=True,
    )


if __name__ == "__main__":
    main()
