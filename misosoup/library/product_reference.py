"""Reference metabolic production profile for microbial communities."""

import logging
from time import perf_counter

from reframed.solvers.solution import Status

from ..reframed.layered_community import BOUND_INF, LayeredCommunity


def get_community_exchanges(community: LayeredCommunity) -> list:
    """Return global community exchange reactions."""
    return [
        rid
        for rid in community.merged_model.reactions
        if rid.startswith("R_EX_") and not rid.endswith("_i")
    ]


def constrain_full_community(
    community: LayeredCommunity,
    minimal_growth: float = 0.01,
) -> None:
    """Force every organism in the community to be active."""
    if not community.has_binary_variables:
        community.setup_binary_variables(minimal_growth)

    y_variables = {
        f"y_{org_id}": 1
        for org_id in community.organisms
    }

    community.solver.add_constraint(
        "c_full_community",
        y_variables,
        "=",
        len(y_variables),
    )
    community.solver.update()


def constrain_full_community_lp(
    community: LayeredCommunity,
    minimal_growth: float = 0.01,
) -> None:
    """Force every organism active with direct bounds and no y_j binaries."""
    if minimal_growth <= 0:
        raise ValueError("minimal_growth must be > 0.")

    if community.has_binary_variables:
        raise ValueError(
            "constrain_full_community_lp requires a community without "
            "organism binary variables."
        )

    fixed_bounds = {}

    for org_id, org_model in community.organisms.items():
        biomass_rid = org_model.biomass_reaction

        for r_id in org_model.reactions:
            if not r_id.startswith("R_EX") and r_id != biomass_rid:
                continue

            merged_id = community.reaction_map[(org_id, r_id)]

            if r_id == biomass_rid:
                fixed_bounds[merged_id] = (minimal_growth, BOUND_INF)
            elif r_id.startswith("R_EX"):
                fixed_bounds[merged_id] = (-BOUND_INF, BOUND_INF)

    gurobi_model = community.solver.problem
    gurobi_model.update()

    missing_variables = []

    for rid, (lower_bound, upper_bound) in fixed_bounds.items():
        variable = gurobi_model.getVarByName(rid)

        if variable is None:
            missing_variables.append(rid)
            continue

        variable.LB = lower_bound
        variable.UB = upper_bound

    if missing_variables:
        raise RuntimeError(
            "Unable to apply full-community LP bounds; missing solver "
            "variables: " + ", ".join(missing_variables[:20])
        )

    gurobi_model.update()

    print(
        "Full-community reference configured as a pure LP: "
        f"{len(community.organisms)} organisms forced active with "
        f"growth >= {minimal_growth:.8g}; "
        f"{len(fixed_bounds)} direct reaction bounds; "
        "0 organism y_j binaries.",
        flush=True,
    )

def find_producible_exchanges(
    community: LayeredCommunity,
    tolerance: float = 1e-6,
    max_exchanges: int | None = None,
    exchange_reactions: list[str] | None = None,
) -> dict:
    """Compute maximum secretion of selected global community exchanges.

    When exchange_reactions is None, every global exchange is considered.
    Callers can pass a biologically prefiltered exchange list to avoid
    expensive optimizations for water, protons, H2, metals, and other
    out-of-scope products.

    max_exchanges is intended for deterministic benchmarks and tests.
    None scans every selected exchange; a positive integer scans only the
    first N exchange IDs after sorting.
    """
    if exchange_reactions is None:
        exchanges = sorted(get_community_exchanges(community))
    else:
        available = set(get_community_exchanges(community))
        requested = set(exchange_reactions)
        unknown = sorted(requested - available)
        if unknown:
            raise ValueError(
                "Requested product-scan exchanges are not global community "
                "exchanges: " + ", ".join(unknown[:20])
            )
        exchanges = sorted(requested)

    if max_exchanges is not None:
        if max_exchanges < 1:
            raise ValueError("max_exchanges must be >= 1 or None.")
        exchanges = exchanges[:max_exchanges]
    producible = {}

    gurobi_model = community.solver.problem
    gurobi_model.update()

    if gurobi_model.NumIntVars != 0:
        raise RuntimeError(
            "Product scan must be a pure LP, but Gurobi reports "
            f"{gurobi_model.NumIntVars} integer variables "
            f"({gurobi_model.NumBinVars} binary)."
        )

    total_exchanges = len(exchanges)
    scan_start = perf_counter()

    if max_exchanges is None:
        scan_label = "full"
    else:
        scan_label = f"benchmark limit={max_exchanges}"

    print(
        f"Product scan: {total_exchanges} exchange reactions selected "
        f"({scan_label}).",
        flush=True,
    )
    print(
        "Product scan LP: "
        f"variables={gurobi_model.NumVars}; "
        f"linear_constraints={gurobi_model.NumConstrs}; "
        f"general_constraints={gurobi_model.NumGenConstrs}; "
        "integer_variables=0; binary_variables=0.",
        flush=True,
    )

    logging.info(
        "Scanning %i community exchange reactions.",
        total_exchanges,
    )

    for index, rid in enumerate(exchanges, start=1):
        solve_start = perf_counter()
        solution = community.solver.solve(
            objective={rid: 1},
            get_values=[rid],
            minimize=False,
        )

        solve_time = perf_counter() - solve_start
        elapsed_time = perf_counter() - scan_start

        if solution.status != Status.OPTIMAL:
            logging.warning(
                "Unable to maximize %s. Solver status: %s",
                rid,
                solution.status,
            )
            print(
                f"Optimization {index}/{total_exchanges} finished: "
                f"{rid} status={solution.status} "
                f"solve={solve_time:.3f}s elapsed={elapsed_time:.3f}s",
                flush=True,
            )
            continue

        maximum = solution.values.get(rid, 0.0)

        if maximum > tolerance:
            producible[rid] = maximum

        print(
            f"Optimization {index}/{total_exchanges} finished: "
            f"{rid} max={maximum:.8g} "
            f"solve={solve_time:.3f}s elapsed={elapsed_time:.3f}s",
            flush=True,
        )

        logging.debug(
            "[%i/%i] %s = %.8g",
            index,
            total_exchanges,
            rid,
            maximum,
        )

    logging.info(
        "Found %i producible metabolites out of %i exchanges.",
        len(producible),
        len(exchanges),
    )

    return producible