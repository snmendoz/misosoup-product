"""Fast one-model blocked-reaction pruning with SciPy/HiGHS.

This implementation is intentionally independent of Gurobi/WLS so hundreds of
models can be processed concurrently with Slurm arrays.

A reaction is blocked iff, with all environmental exchanges open and biomass
bounded to [0, open_bound], it cannot carry flux in either direction. We use an
exact blockedness-oriented FVA: for irreversible reactions only the feasible
direction is optimized; for reversible reactions one direction is optimized
first and the opposite direction is solved only when needed to prove the
reaction blocked. This gives the same blocked/unblocked classification as full
two-sided FVA while avoiding unnecessary LP solves.
"""

from __future__ import annotations

import argparse
import math
import os
from pathlib import Path
from time import perf_counter

import numpy as np
import scipy
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from reframed import ReactionType, save_cbmodel

from common import yaml_dump_atomic, yaml_load
from misosoup.library.readwrite import load_models


SOLVER_BACKEND = "scipy-highs-blockedness-fva-v1"


def source_fingerprint(path: Path) -> dict:
    stat = path.stat()
    return {
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }


def audit_path_for(output_model: Path) -> Path:
    return output_model.with_name(
        output_model.stem + "_blocked_reactions.yaml"
    )


def environmental_exchange_ids(model) -> list[str]:
    """Return environmental exchanges, excluding the biomass drain."""
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


def _finite_or_none(value):
    if value is None:
        return None

    numeric = float(value)
    if math.isinf(numeric):
        return None

    return numeric


def build_lp(model, open_bound: float):
    reaction_ids = list(model.reactions.keys())
    metabolite_ids = list(model.metabolites.keys())

    reaction_index = {
        rid: index
        for index, rid in enumerate(reaction_ids)
    }
    metabolite_index = {
        mid: index
        for index, mid in enumerate(metabolite_ids)
    }

    rows = []
    cols = []
    data = []

    for rid, reaction in model.reactions.items():
        column = reaction_index[rid]

        for mid, coefficient in reaction.stoichiometry.items():
            rows.append(metabolite_index[mid])
            cols.append(column)
            data.append(float(coefficient))

    matrix = coo_matrix(
        (data, (rows, cols)),
        shape=(len(metabolite_ids), len(reaction_ids)),
        dtype=float,
    ).tocsr()

    exchanges = set(environmental_exchange_ids(model))
    biomass = model.biomass_reaction
    bounds = []

    for rid in reaction_ids:
        reaction = model.reactions[rid]

        if rid == biomass:
            bounds.append((0.0, float(open_bound)))
        elif rid in exchanges:
            bounds.append(
                (-float(open_bound), float(open_bound))
            )
        else:
            bounds.append(
                (
                    _finite_or_none(reaction.lb),
                    _finite_or_none(reaction.ub),
                )
            )

    return {
        "reaction_ids": reaction_ids,
        "reaction_index": reaction_index,
        "matrix": matrix,
        "rhs": np.zeros(len(metabolite_ids), dtype=float),
        "bounds": bounds,
        "exchange_ids": sorted(exchanges),
    }


def solve_objective(lp, reaction_index: int, direction: str):
    objective = np.zeros(len(lp["reaction_ids"]), dtype=float)

    if direction == "max":
        objective[reaction_index] = -1.0
    elif direction == "min":
        objective[reaction_index] = 1.0
    else:
        raise ValueError(f"Unknown LP direction: {direction}")

    result = linprog(
        objective,
        A_eq=lp["matrix"],
        b_eq=lp["rhs"],
        bounds=lp["bounds"],
        method="highs",
        options={"presolve": True},
    )

    if not result.success:
        raise RuntimeError(
            f"HiGHS LP failed: status={result.status}; "
            f"message={result.message}"
        )

    flux = float(result.x[reaction_index])
    return flux


def maximum_biomass(model, lp) -> float:
    index = lp["reaction_index"][model.biomass_reaction]
    return solve_objective(lp, index, "max")


def blocked_reactions(
    model,
    lp,
    tolerance: float,
    progress_every: int = 250,
) -> tuple[list[str], int]:
    """Return exact blocked set using the minimum LP solves needed.

    A reaction is blocked iff max(v_r) <= tol AND min(v_r) >= -tol.
    Bounds let us skip one direction for irreversible reactions, and any
    nonzero feasible optimum proves a reaction is unblocked immediately.
    """
    blocked = []
    lp_solves = 0
    total = len(lp["reaction_ids"])
    start = perf_counter()

    for position, rid in enumerate(lp["reaction_ids"], start=1):
        index = lp["reaction_index"][rid]
        lower, upper = lp["bounds"][index]

        lower_numeric = -math.inf if lower is None else float(lower)
        upper_numeric = math.inf if upper is None else float(upper)

        # Fixed zero is trivially blocked.
        if (
            abs(lower_numeric) <= tolerance
            and abs(upper_numeric) <= tolerance
        ):
            blocked.append(rid)
        elif lower_numeric >= -tolerance:
            # Nonnegative reaction: only maximum is needed.
            maximum = solve_objective(lp, index, "max")
            lp_solves += 1

            if maximum <= tolerance:
                blocked.append(rid)
        elif upper_numeric <= tolerance:
            # Nonpositive reaction: only minimum is needed.
            minimum = solve_objective(lp, index, "min")
            lp_solves += 1

            if minimum >= -tolerance:
                blocked.append(rid)
        else:
            # Reversible: a positive maximum proves unblocked. Only when the
            # positive direction cannot carry flux do we test the negative one.
            maximum = solve_objective(lp, index, "max")
            lp_solves += 1

            if maximum <= tolerance:
                minimum = solve_objective(lp, index, "min")
                lp_solves += 1

                if minimum >= -tolerance:
                    blocked.append(rid)

        if (
            progress_every > 0
            and (
                position % progress_every == 0
                or position == total
            )
        ):
            print(
                f"FVA progress: {position}/{total}; "
                f"blocked={len(blocked)}; LP_solves={lp_solves}; "
                f"elapsed={perf_counter() - start:.1f}s",
                flush=True,
            )

    return sorted(blocked), lp_solves


def remove_orphan_metabolites(model) -> list[str]:
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


def cache_valid(
    source_path: Path,
    output_model: Path,
    audit_path: Path,
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
        and audit.get("solver_backend") == SOLVER_BACKEND
        and audit.get("source_model") == str(source_path)
        and audit.get("source_fingerprint")
        == source_fingerprint(source_path)
        and float(audit.get("blocked_tolerance", -1.0))
        == float(tolerance)
        and float(audit.get("open_exchange_bound", -1.0))
        == float(open_bound)
        and audit.get("output_model") == str(output_model)
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-model", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--blocked-tolerance", type=float, default=1e-9)
    parser.add_argument("--open-exchange-bound", type=float, default=1000.0)
    parser.add_argument("--progress-every", type=int, default=250)
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    source_path = args.source_model.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not source_path.exists():
        raise FileNotFoundError(source_path)
    if args.blocked_tolerance <= 0:
        raise ValueError("--blocked-tolerance must be > 0.")
    if args.open_exchange_bound <= 0:
        raise ValueError("--open-exchange-bound must be > 0.")

    output_model = output_dir / (
        source_path.stem + "_without_blocked_reactions.xml"
    )
    audit_path = audit_path_for(output_model)

    if cache_valid(
        source_path,
        output_model,
        audit_path,
        args.blocked_tolerance,
        args.open_exchange_bound,
    ):
        print(
            f"CACHED: {source_path.name} -> {output_model.name}",
            flush=True,
        )
        return

    total_start = perf_counter()
    print(
        f"START {source_path.name}; backend={SOLVER_BACKEND}; "
        f"scipy={scipy.__version__}",
        flush=True,
    )

    model = load_models([str(source_path)])[0]
    biomass = model.biomass_reaction

    if biomass not in model.reactions:
        raise RuntimeError(
            f"Biomass reaction {biomass!r} absent from {source_path}."
        )

    original_reactions = len(model.reactions)
    original_metabolites = len(model.metabolites)

    lp = build_lp(model, args.open_exchange_bound)
    biomass_before = maximum_biomass(model, lp)

    fva_start = perf_counter()
    blocked, lp_solves = blocked_reactions(
        model,
        lp,
        args.blocked_tolerance,
        progress_every=args.progress_every,
    )
    fva_seconds = perf_counter() - fva_start

    if biomass in blocked:
        raise RuntimeError(
            f"Biomass {biomass} classified as blocked with all "
            "environmental exchanges open."
        )

    for rid in blocked:
        model.remove_reaction(rid)

    orphan_metabolites = remove_orphan_metabolites(model)

    reduced_lp = build_lp(model, args.open_exchange_bound)
    biomass_after = maximum_biomass(model, reduced_lp)

    comparison_tolerance = max(
        args.blocked_tolerance,
        1e-7 * max(1.0, abs(biomass_before)),
    )
    if abs(biomass_after - biomass_before) > comparison_tolerance:
        raise RuntimeError(
            "Biomass changed after pruning: "
            f"before={biomass_before:.12g}; "
            f"after={biomass_after:.12g}; "
            f"allowed={comparison_tolerance:.12g}."
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
        "solver_backend": SOLVER_BACKEND,
        "scipy_version": scipy.__version__,
        "source_model": str(source_path),
        "source_fingerprint": source_fingerprint(source_path),
        "output_model": str(output_model),
        "model_id": str(model.id),
        "biomass_reaction": str(biomass),
        "blocked_tolerance": float(args.blocked_tolerance),
        "open_exchange_bound": float(args.open_exchange_bound),
        "biomass_lower_bound_fva": 0.0,
        "biomass_upper_bound_fva": float(
            args.open_exchange_bound
        ),
        "number_environmental_exchanges_opened": len(
            lp["exchange_ids"]
        ),
        "original_reactions": int(original_reactions),
        "blocked_reactions": int(len(blocked)),
        "remaining_reactions": int(len(model.reactions)),
        "reaction_reduction_fraction": (
            len(blocked) / original_reactions
            if original_reactions
            else 0.0
        ),
        "original_metabolites": int(original_metabolites),
        "orphan_metabolites_removed": int(
            len(orphan_metabolites)
        ),
        "remaining_metabolites": int(len(model.metabolites)),
        "biomass_max_open_exchanges_before": float(
            biomass_before
        ),
        "biomass_max_open_exchanges_after": float(
            biomass_after
        ),
        "fva_lp_solves": int(lp_solves),
        "fva_seconds": float(fva_seconds),
        "blocked_reaction_ids": blocked,
        "orphan_metabolite_ids": orphan_metabolites,
        "elapsed_seconds": float(
            perf_counter() - total_start
        ),
    }
    yaml_dump_atomic(audit, audit_path)

    print(
        f"COMPLETE {source_path.name}: "
        f"{original_reactions} -> {len(model.reactions)} reactions "
        f"(-{len(blocked)}, "
        f"{100.0 * len(blocked) / original_reactions:.1f}%); "
        f"LP_solves={lp_solves}; "
        f"FVA={fva_seconds:.1f}s; "
        f"total={perf_counter() - total_start:.1f}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
