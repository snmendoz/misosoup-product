"""Minimum-cardinality communities that preserve metabolic product production.

This module implements Stage D of the product-preservation workflow.

Reference phenotype
===================

Stages A-C define a full-community reference phenotype and a selected product
set S_B.  Stage C provides a parsimonious reference flux:

    v_i^ref

for every selected product i.

Stage D asks:

    What is the smallest microbial community that can preserve those products?

Mathematical formulation
========================

For every organism j introduce the existing MiSoSoup activity variable:

    y_j in {0, 1}

where y_j = 1 means that organism j is active.

For each selected product i, require:

    v_i >= beta * v_i^ref

where:

    beta = product_retention

The minimum-community problem is then:

    min  sum_j y_j

subject to:

    S v = 0

    medium constraints

    y_j = 1  =>  biomass_j >= mu_min

    y_j = 0  =>  organism exchange fluxes are disabled
                  by LayeredCommunity.setup_binary_variables()

    v_i >= beta * v_i^ref       for all i in S_B

The first optimization gives the minimum cardinality:

    N_min = min sum_j y_j

To enumerate all communities with this minimum cardinality, we then impose:

    sum_j y_j = N_min

and repeatedly add a no-good cut for every solution C_k:

    sum_(j in C_k) y_j <= N_min - 1

Because the total cardinality is fixed to N_min, this excludes exactly the
previously found organism combination while retaining all other combinations
of the same size.
"""

import logging
from pathlib import Path

from reframed.solvers.solution import Status

from ..reframed.layered_community import LayeredCommunity


def _safe_name(value: str) -> str:
    """Return a solver-safe name fragment."""

    return "".join(
        character if character.isalnum() or character == "_" else "_"
        for character in value
    )


def _raise_with_iis(
    community: LayeredCommunity,
    stage_name: str,
    filename: str,
    status,
) -> None:
    """Compute and write a Gurobi IIS before raising an error."""

    gurobi_model = community.solver.problem

    print()
    print("=" * 70)
    print(f"{stage_name} IS INFEASIBLE - COMPUTING IIS")
    print("=" * 70)

    try:
        gurobi_model.computeIIS()

        iis_path = Path(filename).resolve()
        gurobi_model.write(str(iis_path))

        print(
            f"IIS written to: {iis_path}",
            flush=True,
        )

    except Exception as error:
        print(
            f"Unable to compute/write IIS: {error}",
            flush=True,
        )

    raise RuntimeError(
        f"Unable to solve {stage_name}. "
        f"Solver status: {status}"
    )


def _selected_organisms(
    solution,
    organism_variables: dict,
    cutoff: float = 0.5,
) -> list:
    """Return organism IDs whose binary variable y_j is active."""

    return [
        org_id
        for org_id, y_name in organism_variables.items()
        if solution.values.get(y_name, 0.0) > cutoff
    ]


def _extract_community_solution(
    community: LayeredCommunity,
    solution,
    organism_variables: dict,
    selected_products: list,
    product_requirements: dict,
    medium: dict,
    tolerance: float,
) -> dict:
    """Convert one minimum-community solution into a serializable dictionary."""

    organisms = _selected_organisms(
        solution,
        organism_variables,
    )

    product_fluxes = {
        rid: float(
            solution.values.get(
                rid,
                0.0,
            )
        )
        for rid in selected_products
    }

    # Validate every preserved product explicitly.
    for rid, required_flux in product_requirements.items():
        observed_flux = product_fluxes[rid]

        if observed_flux + tolerance < required_flux:
            raise RuntimeError(
                f"Minimum-community solution violates product requirement "
                f"for {rid}: observed={observed_flux}, "
                f"required>={required_flux}."
            )

    organism_growth = {}

    for org_id in organisms:
        org_model = community.organisms[org_id]

        biomass_rid = community.reaction_map[
            (
                org_id,
                org_model.biomass_reaction,
            )
        ]

        organism_growth[org_id] = float(
            solution.values.get(
                biomass_rid,
                0.0,
            )
        )

    community_growth_rid = (
        community.merged_model.biomass_reaction
    )

    community_growth = float(
        solution.values.get(
            community_growth_rid,
            0.0,
        )
    )

    # Preserve the raw global-exchange fluxes for every compound in the
    # configured medium.  Negative exchange flux is uptake, while positive
    # exchange flux is secretion.  These values are used by the experiment
    # layer to audit whether a minimum community is relying on uptake bounds.
    medium_exchange_fluxes = {
        rid: float(
            solution.values.get(
                rid,
                0.0,
            )
        )
        for rid in medium
    }

    return {
        "organisms": organisms,
        "size": len(organisms),
        "organism_growth": organism_growth,
        "community_growth": community_growth,
        "product_fluxes": product_fluxes,
        "medium_exchange_fluxes": medium_exchange_fluxes,
    }


def findMinimalProductCommunities(
    community: LayeredCommunity,
    medium: dict,
    reference_products: dict,
    minimal_growth: float = 0.01,
    product_retention: float = 0.90,
    tolerance: float = 1e-6,
    max_communities: int | None = 100,
) -> dict:
    """Find and enumerate minimum-cardinality product-preserving communities.

    Parameters
    ----------
    community
        A FRESH LayeredCommunity containing all candidate organisms.

        It must not contain the ``c_full_community`` constraint from the
        reference calculation, because Stage D needs the y_j variables to
        be free so the optimizer can remove organisms.

    medium
        Community medium dictionary.

    reference_products
        Stage-C ``reference_products`` dictionary.  Each product must contain:

            reference_products[rid]["reference_flux"]

    minimal_growth
        Minimum growth required for every active organism:

            y_j = 1  =>  biomass_j >= minimal_growth

    product_retention
        beta in:

            v_i >= beta * v_i^ref

        Example:

            beta = 0.90

        requires at least 90% of every Stage-C reference product flux.

    tolerance
        Numerical tolerance used when validating returned product fluxes.

    max_communities
        Maximum number of minimum-cardinality communities to enumerate.

        ``None`` means exhaustive enumeration.  A finite cap is useful for
        large communities because the number of alternative minimum
        communities can be combinatorial.

    Returns
    -------
    dict
        {
            "minimum_size": N_min,
            "product_retention": beta,
            "product_requirements": {...},
            "communities": [...],
            "number_communities": ...,
            "enumeration_complete": bool,
        }
    """

    # ==================================================================
    # INPUT VALIDATION
    # ==================================================================

    if not 0 < product_retention <= 1:
        raise ValueError(
            "product_retention must be > 0 and <= 1."
        )

    if minimal_growth <= 0:
        raise ValueError(
            "minimal_growth must be > 0."
        )

    if tolerance < 0:
        raise ValueError(
            "tolerance must be >= 0."
        )

    if max_communities is not None and max_communities < 1:
        raise ValueError(
            "max_communities must be >= 1 or None."
        )

    if not reference_products:
        raise ValueError(
            "reference_products is empty."
        )

    print(
        f"Stage D.0: fresh community with {len(community.organisms)} candidate organisms; "
        f"preserving {len(reference_products)} products at "
        f"{100.0 * product_retention:.3g}% of Stage-C reference flux.",
        flush=True,
    )

    # ==================================================================
    # D.0 CREATE ORGANISM-SELECTION VARIABLES AND APPLY MEDIUM
    # ==================================================================
    #
    # LayeredCommunity.setup_binary_variables() introduces:
    #
    #     y_j in {0,1}
    #
    # and couples organism activity to growth/exchange fluxes.
    #
    # In particular:
    #
    #     y_j = 1  => biomass_j >= minimal_growth
    #
    # while organism exchange reactions are forced off when y_j = 0.
    # ==================================================================

    if not community.has_binary_variables:
        community.setup_binary_variables(
            minimal_growth
        )

    community.setup_medium(
        medium
    )

    organism_variables = {
        org_id: f"y_{org_id}"
        for org_id in community.organisms
    }

    print(
        f"Stage D.0 COMPLETE: {len(organism_variables)} organism-selection "
        "variables and medium constraints installed.",
        flush=True,
    )

    print(
        "Stage D.1: adding product-preservation constraints...",
        flush=True,
    )

    # ==================================================================
    # D.1 PRODUCT-PRESERVATION CONSTRAINTS
    # ==================================================================
    #
    # For each Stage-C selected product:
    #
    #     v_i >= beta * v_i^ref
    #
    # The requirement is defined individually for every product, rather
    # than only preserving a summed score.  Thus no selected product can
    # disappear while another product compensates for it.
    # ==================================================================

    product_requirements = {}

    for index, (rid, info) in enumerate(
        reference_products.items()
    ):
        if rid not in community.merged_model.reactions:
            raise KeyError(
                f"Reference product {rid} is absent from "
                "the minimum-community model."
            )

        reference_flux = float(
            info["reference_flux"]
        )

        if reference_flux <= tolerance:
            raise ValueError(
                f"Reference product {rid} has non-positive/"
                f"negligible reference flux {reference_flux}."
            )

        required_flux = (
            product_retention
            * reference_flux
        )

        product_requirements[rid] = (
            required_flux
        )

        community.solver.add_constraint(
            (
                f"c_mincomm_product_"
                f"{index}_"
                f"{_safe_name(rid)}"
            ),
            {
                rid: 1,
            },
            ">",
            required_flux,
        )

    community.solver.update()

    selected_products = list(
        product_requirements.keys()
    )

    print(
        f"Stage D.1 COMPLETE: {len(selected_products)} product constraints added.",
        flush=True,
    )

    # ==================================================================
    # VALUES TO RETRIEVE FROM EVERY SOLUTION
    # ==================================================================

    biomass_reactions = []

    for org_id, org_model in community.organisms.items():
        biomass_reactions.append(
            community.reaction_map[
                (
                    org_id,
                    org_model.biomass_reaction,
                )
            ]
        )

    # Only request medium exchanges that actually exist in the merged
    # community.  A medium can contain compounds absent from a particular
    # community; setup_medium() warns and skips those reactions.  ReFramed's
    # Gurobi backend cannot retrieve values for non-existent variable names.
    available_medium_reactions = [
        rid
        for rid in medium
        if rid in community.merged_model.reactions
    ]

    values_to_get = (
        list(organism_variables.values())
        + selected_products
        + available_medium_reactions
        + biomass_reactions
        + [community.merged_model.biomass_reaction]
    )

    # ==================================================================
    # D.2 MINIMIZE COMMUNITY CARDINALITY
    # ==================================================================
    #
    # Objective:
    #
    #     N_min = min sum_j y_j
    # ==================================================================

    community_objective = {
        y_name: 1
        for y_name in organism_variables.values()
    }

    logging.info(
        "Stage D: minimizing community size while preserving %i products.",
        len(selected_products),
    )

    print(
        f"Stage D.2: solving minimum-cardinality MILP for "
        f"{len(organism_variables)} organisms and {len(selected_products)} products...",
        flush=True,
    )

    first_solution = community.solver.solve(
        objective=community_objective,
        get_values=values_to_get,
        minimize=True,
    )

    if first_solution.status != Status.OPTIMAL:
        _raise_with_iis(
            community=community,
            stage_name="Stage D minimum-community optimization",
            filename="stage_d_minimum_community_iis.ilp",
            status=first_solution.status,
        )

    first_selected = _selected_organisms(
        first_solution,
        organism_variables,
    )

    minimum_size = len(
        first_selected
    )

    if minimum_size == 0:
        raise RuntimeError(
            "Stage D returned a zero-organism community despite "
            "positive product requirements."
        )

    logging.info(
        "Stage D minimum community size N_min = %i.",
        minimum_size,
    )

    print(
        f"Stage D.2 COMPLETE: global minimum community size N_min={minimum_size}.",
        flush=True,
    )

    # ==================================================================
    # D.3 FIX CARDINALITY TO N_min
    # ==================================================================
    #
    # We now enumerate ONLY globally minimum-cardinality communities:
    #
    #     sum_j y_j = N_min
    # ==================================================================

    community.solver.add_constraint(
        "c_mincomm_fix_minimum_size",
        community_objective,
        "=",
        minimum_size,
    )

    community.solver.update()

    print(
        f"Stage D.3: cardinality fixed at N_min={minimum_size}; "
        "recording first optimum and preparing alternative enumeration.",
        flush=True,
    )

    communities = [
        _extract_community_solution(
            community=community,
            solution=first_solution,
            organism_variables=organism_variables,
            selected_products=selected_products,
            product_requirements=product_requirements,
            medium=medium,
            tolerance=tolerance,
        )
    ]

    # ==================================================================
    # D.4 ENUMERATE ALTERNATIVE MINIMUM COMMUNITIES
    # ==================================================================
    #
    # For a previously found community C_k with |C_k| = N_min:
    #
    #     sum_(j in C_k) y_j <= N_min - 1
    #
    # Together with:
    #
    #     sum_j y_j = N_min
    #
    # this excludes exactly that organism combination.
    # ==================================================================

    no_good_index = 0

    def add_no_good_cut(selected_organisms: list) -> None:
        nonlocal no_good_index

        community.solver.add_constraint(
            f"c_mincomm_nogood_{no_good_index}",
            {
                organism_variables[org_id]: 1
                for org_id in selected_organisms
            },
            "<",
            minimum_size - 1,
        )

        no_good_index += 1
        community.solver.update()

    add_no_good_cut(
        first_selected
    )

    enumeration_complete = False

    print(
        f"Stage D.4: enumerating alternative minimum communities "
        f"(cap={max_communities if max_communities is not None else 'none'}).",
        flush=True,
    )
    print(
        f"Stage D.4: community 1 found; size={minimum_size}; members={first_selected}.",
        flush=True,
    )

    while True:

        # Stop early if the user requested an enumeration cap.
        if (
            max_communities is not None
            and len(communities) >= max_communities
        ):
            logging.warning(
                "Stage D reached max_communities=%i before "
                "proving enumeration complete.",
                max_communities,
            )
            break

        # With cardinality already fixed, a zero objective is sufficient.
        print(
            f"Stage D.4: searching for community {len(communities) + 1}...",
            flush=True,
        )
        solution = community.solver.solve(
            objective={},
            get_values=values_to_get,
            minimize=True,
        )

        # No more feasible solution of size N_min:
        # enumeration is complete only when the solver explicitly proves
        # infeasibility. Other statuses (for example numerical failure or a
        # time limit) must not be silently misclassified as exhaustive
        # enumeration.
        if solution.status == Status.INFEASIBLE:
            enumeration_complete = True
            break

        if solution.status != Status.OPTIMAL:
            raise RuntimeError(
                "Stage D enumeration stopped before completion. "
                f"Solver status: {solution.status}"
            )

        selected = _selected_organisms(
            solution,
            organism_variables,
        )

        if len(selected) != minimum_size:
            raise RuntimeError(
                "Stage D enumeration returned a community with "
                f"size {len(selected)} although N_min={minimum_size}."
            )

        communities.append(
            _extract_community_solution(
                community=community,
                solution=solution,
                organism_variables=organism_variables,
                selected_products=selected_products,
                product_requirements=product_requirements,
                tolerance=tolerance,
            )
        )

        logging.info(
            "Stage D found minimum community %i: %s",
            len(communities),
            selected,
        )

        print(
            f"Stage D.4: community {len(communities)} found; "
            f"size={minimum_size}; members={selected}.",
            flush=True,
        )

        add_no_good_cut(
            selected
        )

    print(
        f"Stage D COMPLETE: N_min={minimum_size}; "
        f"communities_found={len(communities)}; "
        f"enumeration_complete={enumeration_complete}.",
        flush=True,
    )

    return {
        "minimum_size": minimum_size,
        "product_retention": product_retention,
        "product_requirements": product_requirements,
        "number_communities": len(communities),
        "enumeration_complete": enumeration_complete,
        "communities": communities,
    }
