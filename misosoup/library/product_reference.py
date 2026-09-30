"""Reference metabolic production profile for microbial communities."""

import logging

from reframed.solvers.solution import Status

from ..reframed.layered_community import LayeredCommunity


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


def find_producible_exchanges(
    community: LayeredCommunity,
    tolerance: float = 1e-6,
) -> dict:
    """Compute maximum secretion of every global community exchange."""
    exchanges = get_community_exchanges(community)
    producible = {}

    logging.info(
        "Scanning %i community exchange reactions.",
        len(exchanges),
    )

    for index, rid in enumerate(exchanges, start=1):
        solution = community.solver.solve(
            objective={rid: 1},
            get_values=[rid],
            minimize=False,
        )

        if solution.status != Status.OPTIMAL:
            logging.warning(
                "Unable to maximize %s. Solver status: %s",
                rid,
                solution.status,
            )
            continue

        maximum = solution.values.get(rid, 0.0)

        if maximum > tolerance:
            producible[rid] = maximum

        logging.debug(
            "[%i/%i] %s = %.8g",
            index,
            len(exchanges),
            rid,
            maximum,
        )

    logging.info(
        "Found %i producible metabolites out of %i exchanges.",
        len(producible),
        len(exchanges),
    )

    return producible