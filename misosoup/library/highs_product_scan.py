"""Parallel pure-LP product scan using SciPy/HiGHS.

This module intentionally avoids Gurobi.  A sample community is merged once,
converted to a sparse steady-state LP once, then individual exchange objectives
are solved concurrently in forked worker processes.  SLURM provides the outer
sample-level parallelism; this module provides the inner exchange-level
parallelism.
"""

from __future__ import annotations

import math
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from time import perf_counter
from typing import Callable

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csc_matrix

from ..cobra.layered_community import BOUND_INF, LayeredCommunity


@dataclass
class ProductScanLP:
    reaction_ids: list[str]
    reaction_index: dict[str, int]
    a_eq: csc_matrix
    b_eq: np.ndarray
    bounds: list[tuple[float | None, float | None]]


_WORKER_LP: ProductScanLP | None = None
_WORKER_TOLERANCE = 1e-6


def _finite_bound(value, *, lower: bool):
    if value is None:
        return None
    value = float(value)
    if math.isinf(value):
        return None
    return value


def build_product_scan_lp(
    community: LayeredCommunity,
    *,
    minimal_growth: float,
    medium: dict[str, float],
) -> ProductScanLP:
    """Build the exact pure-LP reference problem used by the Product Scan."""
    if community.solver is not None:
        raise ValueError(
            "build_product_scan_lp expects a merge-only LayeredCommunity "
            "(create_solver=False)."
        )
    if minimal_growth <= 0:
        raise ValueError("minimal_growth must be > 0.")

    model = community.merged_model
    reaction_ids = list(model.reactions.keys())
    metabolite_ids = list(model.metabolites.keys())
    reaction_index = {rid: i for i, rid in enumerate(reaction_ids)}
    metabolite_index = {mid: i for i, mid in enumerate(metabolite_ids)}

    rows = []
    cols = []
    data = []

    lower = np.empty(len(reaction_ids), dtype=float)
    upper = np.empty(len(reaction_ids), dtype=float)

    for j, rid in enumerate(reaction_ids):
        reaction = model.reactions[rid]
        lb = float(reaction.lb)
        ub = float(reaction.ub)
        lower[j] = lb
        upper[j] = ub

        for mid, coefficient in reaction.stoichiometry.items():
            rows.append(metabolite_index[mid])
            cols.append(j)
            data.append(float(coefficient))

    # Match constrain_full_community_lp exactly: every organism is active;
    # biomass has a positive lower bound and local exchanges are bounded.
    for org_id, org_model in community.organisms.items():
        biomass_rid = org_model.biomass_reaction

        for rid in org_model.reactions:
            if not rid.startswith("R_EX") and rid != biomass_rid:
                continue

            merged_rid = community.reaction_map[(org_id, rid)]
            j = reaction_index[merged_rid]

            if rid == biomass_rid:
                lower[j] = float(minimal_growth)
                upper[j] = float(BOUND_INF)
            else:
                lower[j] = -float(BOUND_INF)
                upper[j] = float(BOUND_INF)

    # Match setup_medium: global exchanges not in the medium cannot be taken up;
    # listed medium exchanges receive their specified negative lower bound.
    for rid in reaction_ids:
        if rid.startswith("R_EX_") and not rid.endswith("_i"):
            j = reaction_index[rid]
            medium_lb = float(medium.get(rid, 0.0))
            lower[j] = max(lower[j], medium_lb)

    invalid = [
        reaction_ids[j]
        for j in range(len(reaction_ids))
        if lower[j] > upper[j]
    ]
    if invalid:
        raise ValueError(
            "Product-scan LP has lower bound > upper bound for: "
            + ", ".join(invalid[:20])
        )

    bounds = [
        (
            _finite_bound(lower[j], lower=True),
            _finite_bound(upper[j], lower=False),
        )
        for j in range(len(reaction_ids))
    ]

    a_eq = csc_matrix(
        (data, (rows, cols)),
        shape=(len(metabolite_ids), len(reaction_ids)),
        dtype=float,
    )
    b_eq = np.zeros(len(metabolite_ids), dtype=float)

    return ProductScanLP(
        reaction_ids=reaction_ids,
        reaction_index=reaction_index,
        a_eq=a_eq,
        b_eq=b_eq,
        bounds=bounds,
    )


def check_product_scan_feasibility(
    lp: ProductScanLP,
    *,
    tolerance: float = 1e-6,
) -> dict:
    """Solve one zero-objective LP before dispatching exchange objectives."""
    result = linprog(
        np.zeros(len(lp.reaction_ids), dtype=float),
        A_eq=lp.a_eq,
        b_eq=lp.b_eq,
        bounds=lp.bounds,
        method="highs",
        options={
            "primal_feasibility_tolerance": float(tolerance),
            "dual_feasibility_tolerance": float(tolerance),
            "presolve": True,
        },
    )
    return {
        "optimal": bool(result.status == 0),
        "status": int(result.status),
        "message": str(result.message),
    }


def _solve_one_exchange(rid: str) -> tuple[str, dict]:
    lp = _WORKER_LP
    if lp is None:
        raise RuntimeError("Product-scan worker LP was not initialized.")

    j = lp.reaction_index[rid]
    objective = np.zeros(len(lp.reaction_ids), dtype=float)
    objective[j] = -1.0

    start = perf_counter()
    result = linprog(
        objective,
        A_eq=lp.a_eq,
        b_eq=lp.b_eq,
        bounds=lp.bounds,
        method="highs",
        options={
            "primal_feasibility_tolerance": float(_WORKER_TOLERANCE),
            "dual_feasibility_tolerance": float(_WORKER_TOLERANCE),
            "presolve": True,
        },
    )
    solve_seconds = perf_counter() - start

    if result.status != 0:
        return rid, {
            "status": "nonoptimal",
            "solver_status": f"HiGHS status={result.status}: {result.message}",
            "maximum": None,
            "solve_seconds": float(solve_seconds),
        }

    maximum = float(result.x[j])
    return rid, {
        "status": "optimal",
        "solver_status": "HiGHS optimal",
        "maximum": maximum,
        "producible": bool(maximum > _WORKER_TOLERANCE),
        "solve_seconds": float(solve_seconds),
    }


def run_parallel_product_scan(
    lp: ProductScanLP,
    exchange_reactions: list[str],
    *,
    workers: int,
    tolerance: float = 1e-6,
    completed_results: dict[str, dict] | None = None,
    progress_callback: Callable[[str, dict], None] | None = None,
) -> tuple[dict[str, float], dict[str, dict], float]:
    """Maximize selected exchanges concurrently with forked HiGHS workers."""
    if workers < 1:
        raise ValueError("workers must be >= 1.")

    requested = sorted(exchange_reactions)
    unknown = sorted(set(requested) - set(lp.reaction_index))
    if unknown:
        raise ValueError(
            "Unknown product-scan exchange reactions: "
            + ", ".join(unknown[:20])
        )

    results = dict(completed_results or {})
    reusable = {
        rid: record
        for rid, record in results.items()
        if rid in requested
        and record.get("status") == "optimal"
        and record.get("maximum") is not None
    }
    remaining = [rid for rid in requested if rid not in reusable]

    producible = {
        rid: float(record["maximum"])
        for rid, record in reusable.items()
        if float(record["maximum"]) > tolerance
    }

    print(
        f"Product Scan resume: {len(reusable)}/{len(requested)} restored; "
        f"{len(remaining)} remaining.",
        flush=True,
    )
    if not remaining:
        return producible, results, 0.0

    global _WORKER_LP, _WORKER_TOLERANCE
    _WORKER_LP = lp
    _WORKER_TOLERANCE = float(tolerance)

    wall_start = perf_counter()

    def accept(rid: str, record: dict) -> None:
        record = dict(record)
        record["elapsed_seconds_current_run"] = float(
            perf_counter() - wall_start
        )
        results[rid] = record
        if (
            record.get("status") == "optimal"
            and record.get("maximum") is not None
            and float(record["maximum"]) > tolerance
        ):
            producible[rid] = float(record["maximum"])
        if progress_callback is not None:
            progress_callback(rid, record)

    if workers == 1:
        for position, rid in enumerate(remaining, start=1):
            print(
                f"Product Scan {position}/{len(remaining)}: {rid}",
                flush=True,
            )
            solved_rid, record = _solve_one_exchange(rid)
            accept(solved_rid, record)
    else:
        if "fork" not in mp.get_all_start_methods():
            raise RuntimeError(
                "Parallel Product Scan requires Linux/fork multiprocessing."
            )

        context = mp.get_context("fork")
        print(
            f"Product Scan: dispatching {len(remaining)} exchanges to "
            f"{workers} concurrent HiGHS worker processes.",
            flush=True,
        )

        with ProcessPoolExecutor(
            max_workers=workers,
            mp_context=context,
        ) as executor:
            futures = {
                executor.submit(_solve_one_exchange, rid): rid
                for rid in remaining
            }
            completed = 0
            for future in as_completed(futures):
                rid, record = future.result()
                completed += 1
                accept(rid, record)
                print(
                    f"Product Scan completed {completed}/{len(remaining)}: "
                    f"{rid} status={record['status']} "
                    f"max={record.get('maximum')} "
                    f"solve={record.get('solve_seconds', 0.0):.3f}s",
                    flush=True,
                )

    return producible, results, float(perf_counter() - wall_start)
