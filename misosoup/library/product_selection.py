"""Optimization of metabolic product breadth, production, and pFBA reference fluxes.

Workflow
========

Stage A
-------
getMaxProduct()

    Determine the maximum number K* of products that can be produced
    simultaneously above a fraction alpha of their individual maximum M_i.

Stage B.1
---------
maximizeProducts()

    Among solutions containing exactly K* active products, maximize the
    total normalized production.

Stage B.2
---------
getSelectedProductsFromProductMaximization()

    Extract and validate the product set selected by Stage B.

Stage C
-------
solvepFBAUsingFixProducts()

    Fix the product set selected by Stage B and obtain a parsimonious
    reference flux distribution using an absolute-value epigraph pFBA
    formulation:

        a_r >= v_r
        a_r >= -v_r
        a_r >= 0

        min sum_r a_r

    At the optimum:

        a_r = |v_r|

    This formulation leaves the feasible region of the biological flux
    variables unchanged while still minimizing total absolute flux.

Important
---------
Stages A and B mutate one reference-community solver and must be called
sequentially on that same ``LayeredCommunity`` object.

Stage C is intentionally different: it must receive a NEW
``LayeredCommunity`` containing the same organisms. Only the numerical
Stage-B result is transferred into Stage C. This isolates the pFBA
reference problem from the Stage-A/Stage-B MILP machinery.
"""

import logging
import re
from pathlib import Path

from gurobipy import GRB
from reframed.solvers.solution import Status
from reframed.solvers.solver import VarType

from ..reframed.layered_community import LayeredCommunity
from .product_reference import constrain_full_community_lp


# ======================================================================
# GENERAL HELPERS
# ======================================================================


def _safe_name(value: str) -> str:
    """Convert a reaction ID into a safe fragment for constraint names."""

    return re.sub(r"[^A-Za-z0-9_]+", "_", value)


def _raise_with_iis(
    community: LayeredCommunity,
    stage_name: str,
    filename: str,
    status,
) -> None:
    """Compute an IIS when Gurobi reports an infeasible optimization.

    IIS = Irreducible Inconsistent Subsystem.

    It identifies a subset of constraints/bounds that is itself
    inconsistent and is useful for diagnosing infeasibility.
    """

    print()
    print("=" * 70)
    print(f"{stage_name} IS INFEASIBLE - COMPUTING IIS")
    print("=" * 70)

    gurobi_model = community.solver.problem

    try:
        gurobi_model.computeIIS()

        linear_constraints = [
            constraint.ConstrName
            for constraint in gurobi_model.getConstrs()
            if constraint.IISConstr
        ]

        general_constraints = [
            constraint.GenConstrName
            for constraint in gurobi_model.getGenConstrs()
            if constraint.IISGenConstr
        ]

        lower_bounds = [
            variable.VarName
            for variable in gurobi_model.getVars()
            if variable.IISLB
        ]

        upper_bounds = [
            variable.VarName
            for variable in gurobi_model.getVars()
            if variable.IISUB
        ]

        print()
        print("Conflicting linear constraints:")
        if linear_constraints:
            for name in linear_constraints:
                print(f"  {name}")
        else:
            print("  None reported.")

        print()
        print("Conflicting indicator/general constraints:")
        if general_constraints:
            for name in general_constraints:
                print(f"  {name}")
        else:
            print("  None reported.")

        print()
        print("Variables with conflicting lower bounds:")
        if lower_bounds:
            for name in lower_bounds:
                print(f"  {name}")
        else:
            print("  None reported.")

        print()
        print("Variables with conflicting upper bounds:")
        if upper_bounds:
            for name in upper_bounds:
                print(f"  {name}")
        else:
            print("  None reported.")

        iis_path = Path(filename).resolve()
        gurobi_model.write(str(iis_path))

        print()
        print(f"IIS written to: {iis_path}")

    except Exception as error:
        print()
        print(f"Unable to compute/write IIS: {error}")

    raise RuntimeError(
        f"Unable to solve {stage_name}. Solver status: {status}"
    )


# ======================================================================
# STAGE A
# ======================================================================


def getMaxProduct(
    community: LayeredCommunity,
    max_secretion: dict,
    fraction: float = 0.20,
    tolerance: float = 1e-6,
    integer_tolerance: float = 1e-9,
) -> dict:
    """Stage A: determine the maximum number of simultaneous products.

    For each candidate product i:

        M_i = maximum secretion of product i when optimized individually.

    Introduce:

        z_i in {0, 1}

    where z_i = 1 means that product i is counted as simultaneously
    produced.

    The production requirement is:

        z_i = 1  =>  v_i >= alpha * M_i

    where:

        alpha = fraction

    This implication is implemented using a native Gurobi indicator
    constraint; no manually chosen big-M constant is used.

    Stage A solves:

        K* = max sum_i z_i

    subject to:

        S v = 0
        medium constraints
        organism growth constraints
        z_i = 1  =>  v_i >= alpha M_i
        z_i in {0,1}

    Returns a state dictionary used by Stage B.
    """

    # Validate input parameters.
    if not 0 < fraction <= 1:
        raise ValueError("fraction must be > 0 and <= 1.")

    if tolerance < 0:
        raise ValueError("tolerance must be >= 0.")

    if not 1e-9 <= integer_tolerance <= 1e-1:
        raise ValueError(
            "integer_tolerance must be between 1e-9 and 1e-1."
        )

    # Candidate product set:
    #
    #     L = {i : M_i > tolerance}
    products = {
        rid: float(maximum)
        for rid, maximum in max_secretion.items()
        if maximum > tolerance
    }

    if not products:
        raise ValueError(
            "No producible products remain after applying tolerance."
        )

    product_ids = sorted(products)

    print(
        f"Stage A.0: {len(product_ids)} candidate products; "
        f"alpha={fraction:.6g}; preparing binary product-selection MILP.",
        flush=True,
    )

    # Access the Gurobi model underlying the ReFramed solver.
    gurobi_model = community.solver.problem

    # Use a strict binary-integrality tolerance.
    gurobi_model.Params.IntFeasTol = integer_tolerance

    logging.info(
        "Stage A: Gurobi IntFeasTol = %.3g",
        integer_tolerance,
    )

    # Create one binary variable z_i for every candidate product.
    product_variables = {}
    thresholds = {}

    for index, rid in enumerate(product_ids):
        z_name = f"z_product_{index}"

        community.solver.add_variable(
            z_name,
            0,
            1,
            vartype=VarType.BINARY,
        )

        product_variables[rid] = z_name

        # T_i = alpha * M_i
        thresholds[rid] = fraction * products[rid]

    community.solver.update()
    gurobi_model.update()

    print(
        f"Stage A.1: created {len(product_variables)} binary z_i variables; "
        "adding product-threshold indicator constraints.",
        flush=True,
    )

    # Native indicator constraints:
    #
    #     z_i = 1  =>  v_i >= alpha * M_i
    #
    # No user-defined big-M is present.
    for index, rid in enumerate(product_ids):
        z_name = product_variables[rid]

        z_var = gurobi_model.getVarByName(z_name)
        flux_var = gurobi_model.getVarByName(rid)

        if z_var is None:
            raise RuntimeError(f"Unable to find variable {z_name}.")

        if flux_var is None:
            raise RuntimeError(
                f"Unable to find exchange reaction {rid}."
            )

        gurobi_model.addGenConstrIndicator(
            z_var,
            1,
            flux_var,
            GRB.GREATER_EQUAL,
            thresholds[rid],
            name=(
                f"ind_stageA_{index}_{_safe_name(rid)}"
            ),
        )

    gurobi_model.update()

    # Objective:
    #
    #     max sum_i z_i
    objective = {
        z_name: 1
        for z_name in product_variables.values()
    }

    values_to_get = (
        list(product_variables.values())
        + product_ids
    )

    logging.info(
        "Stage A: solving max sum(z_i) for %i products.",
        len(product_ids),
    )

    print(
        f"Stage A.2: solving max sum(z_i) for {len(product_ids)} products...",
        flush=True,
    )

    solution = community.solver.solve(
        objective=objective,
        get_values=values_to_get,
        minimize=False,
    )

    if solution.status != Status.OPTIMAL:
        _raise_with_iis(
            community=community,
            stage_name="Stage A",
            filename="stage_a_iis.ilp",
            status=solution.status,
        )

    # Extract products with z_i = 1.
    selected_products = [
        rid
        for rid in product_ids
        if solution.values.get(
            product_variables[rid],
            0.0,
        ) > 0.5
    ]

    k_star = len(selected_products)

    logging.info(
        "Stage A optimum K*: %i / %i",
        k_star,
        len(product_ids),
    )

    print(
        f"Stage A COMPLETE: K*={k_star}/{len(product_ids)} simultaneous products.",
        flush=True,
    )

    return {
        "products": products,
        "product_ids": product_ids,
        "thresholds": thresholds,
        "product_variables": product_variables,
        "fraction": fraction,
        "tolerance": tolerance,
        "integer_tolerance": integer_tolerance,
        "k_star": k_star,
        "stage_a_selected_products": selected_products,
        "stage_a_solution": solution,
    }


# ======================================================================
# STAGE B.1
# ======================================================================


def maximizeProducts(
    community: LayeredCommunity,
    stage_a: dict,
) -> dict:
    """Stage B.1: maximize normalized production at fixed K*.

    Stage A determined:

        K* = maximum number of simultaneous products.

    Stage B fixes only the NUMBER of active products:

        sum_i z_i = K*

    but still allows Gurobi to choose WHICH K* products are active.

    For each product introduce:

        0 <= q_i <= 1

    with:

        q_i <= z_i

    and:

        z_i = 1  =>  v_i >= M_i q_i

    Therefore, for an active product:

        q_i <= v_i / M_i

    Stage B solves:

        Q* = max sum_i q_i

    subject to:

        sum_i z_i = K*
        q_i <= z_i
        z_i = 1 => v_i >= M_i q_i

    plus all Stage A and metabolic constraints.

    Returns a state dictionary used by Stage B.2 and Stage C.
    """

    products = stage_a["products"]
    product_ids = stage_a["product_ids"]
    product_variables = stage_a["product_variables"]
    k_star = stage_a["k_star"]

    gurobi_model = community.solver.problem

    print(
        f"Stage B.0: fixing product cardinality at K*={k_star}; "
        f"preparing normalized-production optimization for {len(product_ids)} products.",
        flush=True,
    )

    # Fix only the cardinality:
    #
    #     sum_i z_i = K*
    #
    # Product identities remain free.
    community.solver.add_constraint(
        "c_stageA_fix_Kstar",
        {
            z_name: 1
            for z_name in product_variables.values()
        },
        "=",
        k_star,
    )

    community.solver.update()

    # Create normalized-production variables:
    #
    #     0 <= q_i <= 1
    normalized_variables = {}

    for index, rid in enumerate(product_ids):
        q_name = f"q_product_{index}"

        community.solver.add_variable(
            q_name,
            0,
            1,
            vartype=VarType.CONTINUOUS,
        )

        normalized_variables[rid] = q_name

    community.solver.update()
    gurobi_model.update()

    # Require:
    #
    #     q_i <= z_i
    #
    # Thus z_i = 0 forces q_i = 0.
    for index, rid in enumerate(product_ids):
        q_name = normalized_variables[rid]
        z_name = product_variables[rid]

        community.solver.add_constraint(
            f"c_stageB_q_active_{index}_{_safe_name(rid)}",
            {
                q_name: 1,
                z_name: -1,
            },
            "<",
            0,
        )

    community.solver.update()
    gurobi_model.update()

    # Native indicator:
    #
    #     z_i = 1  =>  v_i >= M_i q_i
    #
    # equivalently, when z_i = 1:
    #
    #     v_i - M_i q_i >= 0
    for index, rid in enumerate(product_ids):
        z_var = gurobi_model.getVarByName(
            product_variables[rid]
        )

        q_var = gurobi_model.getVarByName(
            normalized_variables[rid]
        )

        flux_var = gurobi_model.getVarByName(rid)

        if z_var is None or q_var is None or flux_var is None:
            raise RuntimeError(
                f"Unable to resolve Stage B variables for {rid}."
            )

        maximum = products[rid]

        expression = (
            flux_var
            - maximum * q_var
        )

        gurobi_model.addGenConstrIndicator(
            z_var,
            1,
            expression,
            GRB.GREATER_EQUAL,
            0.0,
            name=(
                f"ind_stageB_{index}_{_safe_name(rid)}"
            ),
        )

    gurobi_model.update()

    # Objective:
    #
    #     max sum_i q_i
    objective = {
        q_name: 1
        for q_name in normalized_variables.values()
    }

    values_to_get = (
        list(product_variables.values())
        + list(normalized_variables.values())
        + product_ids
    )

    logging.info(
        "Stage B: maximizing normalized production for K* = %i.",
        k_star,
    )

    print(
        f"Stage B.1: constraints ready; solving max sum(q_i) with K*={k_star}...",
        flush=True,
    )

    solution = community.solver.solve(
        objective=objective,
        get_values=values_to_get,
        minimize=False,
    )

    if solution.status != Status.OPTIMAL:
        _raise_with_iis(
            community=community,
            stage_name="Stage B",
            filename="stage_b_iis.ilp",
            status=solution.status,
        )

    # Q* directly from the solver objective.
    q_star = float(solution.fobj)

    # Independent numerical consistency check.
    q_star_from_values = sum(
        solution.values.get(
            q_name,
            0.0,
        )
        for q_name in normalized_variables.values()
    )

    logging.info(
        "Stage B optimum Q*: %.12g",
        q_star,
    )

    logging.info(
        "Stage B Q* reconstructed from q_i: %.12g",
        q_star_from_values,
    )

    print(
        f"Stage B.1 COMPLETE: Q*={q_star:.12g}; "
        f"reconstructed={q_star_from_values:.12g}.",
        flush=True,
    )

    result = dict(stage_a)

    result.update(
        {
            "normalized_variables": normalized_variables,
            "q_star": q_star,
            "q_star_from_values": q_star_from_values,
            "q_star_difference": (
                q_star
                - q_star_from_values
            ),
            "stage_b_solution": solution,
        }
    )

    return result


# ======================================================================
# STAGE B.2
# ======================================================================


def getSelectedProductsFromProductMaximization(
    stage_b: dict,
    binary_cutoff: float = 0.5,
) -> dict:
    """Stage B.2: extract the product set chosen in Stage B.

    This function performs no optimization.

    It defines:

        S_B = {i : z_i = 1}

    and records, for every selected product:

        M_i
        alpha M_i
        Stage-B flux v_i
        Stage-B q_i
        Stage-B normalized flux v_i/M_i

    The cardinality is validated against:

        |S_B| = K*
    """

    products = stage_b["products"]
    product_ids = stage_b["product_ids"]
    thresholds = stage_b["thresholds"]
    product_variables = stage_b["product_variables"]
    normalized_variables = stage_b["normalized_variables"]
    solution = stage_b["stage_b_solution"]
    k_star = stage_b["k_star"]

    print(
        "Stage B.2: extracting and validating the selected product set...",
        flush=True,
    )

    selected_products = []
    selected_product_info = {}

    for rid in product_ids:
        z_value = solution.values.get(
            product_variables[rid],
            0.0,
        )

        if z_value <= binary_cutoff:
            continue

        selected_products.append(rid)

        flux = solution.values.get(
            rid,
            0.0,
        )

        q_value = solution.values.get(
            normalized_variables[rid],
            0.0,
        )

        selected_product_info[rid] = {
            "max_individual": products[rid],
            "required_flux": thresholds[rid],
            "stage_b_flux": flux,
            "stage_b_q": q_value,
            "stage_b_normalized_flux": (
                flux / products[rid]
            ),
        }

    if len(selected_products) != k_star:
        raise RuntimeError(
            "Stage B product extraction is inconsistent: "
            f"found {len(selected_products)} selected products "
            f"but K* = {k_star}."
        )

    logging.info(
        "Stage B.2: extracted %i selected products.",
        len(selected_products),
    )

    print(
        f"Stage B COMPLETE: selected {len(selected_products)} products "
        f"(expected K*={k_star}).",
        flush=True,
    )

    return {
        "selected_products": selected_products,
        "selected_product_info": selected_product_info,
    }


# ======================================================================
# STAGE C
# ======================================================================


def solvepFBAUsingFixProducts(
    community: LayeredCommunity,
    stage_b: dict,
    product_selection: dict,
    medium: dict,
    minimal_growth: float = 0.01,
    lexicographic_tolerance: float = 1e-5,
    tolerance: float = 1e-6,
    check_feasibility: bool = True,
) -> dict:
    """Stage C: rebuild the selected-product reference problem and run pFBA.

    This function deliberately receives a FRESH ``LayeredCommunity``.

    Stages A/B determine the biological reference quantities:

        S_B     selected product identities
        M_i     individual maximum secretion
        alpha   Stage-A production fraction
        Q*      optimum of the Stage-B normalized-production objective

    Stage C transfers only those numerical results. It does NOT inherit
    product z_i binaries, Stage-A indicators, Stage-B q_i variables,
    Stage-B indicators, or previous Stage-C constraints.

    Diagnostic checkpoints:

        C0: fresh biological problem before pFBA auxiliaries
        C1: absolute-value epigraph auxiliaries added, no pFBA objective
        C2: pFBA minimization

    The formulation is:

        biomass_j >= minimal_growth for every organism j
        no organism-selection binaries
        v_i >= alpha M_i
        0 <= q_i <= 1
        v_i >= M_i q_i
        sum_i q_i >= Q* - epsilon

        a_r >= v_r
        a_r >= -v_r
        a_r >= 0

        min sum_r a_r
    """

    if lexicographic_tolerance < 0:
        raise ValueError(
            "lexicographic_tolerance must be >= 0."
        )

    if minimal_growth <= 0:
        raise ValueError(
            "minimal_growth must be > 0."
        )

    if tolerance < 0:
        raise ValueError(
            "tolerance must be >= 0."
        )

    if community.has_binary_variables:
        raise ValueError(
            "Stage C requires a fresh LayeredCommunity. "
            "The supplied community already has binary variables."
        )

    products = stage_b["products"]
    thresholds = stage_b["thresholds"]
    q_star = float(stage_b["q_star"])
    k_star = int(stage_b["k_star"])

    selected_products = list(
        product_selection["selected_products"]
    )

    if len(selected_products) != k_star:
        raise RuntimeError(
            "Stage C received an unexpected product count: "
            f"{len(selected_products)} instead of K*={k_star}."
        )

    if len(set(selected_products)) != len(selected_products):
        raise RuntimeError(
            "Stage C received duplicated product IDs."
        )

    print(
        f"Stage C.0: fresh community with {len(community.organisms)} organisms; "
        f"{len(selected_products)} selected products; preparing pure-LP pFBA "
        "reference problem.",
        flush=True,
    )

    # Stage C uses the full reference community. Since every organism is
    # fixed active a priori, no organism-selection binaries are needed.
    # Apply the y_j=1 state directly as continuous reaction bounds.
    constrain_full_community_lp(
        community,
        minimal_growth=minimal_growth,
    )

    community.setup_medium(
        medium
    )

    print(
        "Stage C.0: full-community activity and medium constraints installed.",
        flush=True,
    )

    # Fresh q_i variables only for products selected by Stage B.
    normalized_variables = {}

    for index, rid in enumerate(selected_products):

        if rid not in community.merged_model.reactions:
            raise KeyError(
                f"Selected product {rid} is absent from the fresh "
                "Stage-C community."
            )

        if rid not in products:
            raise KeyError(
                f"Selected product {rid} is absent from Stage-B products."
            )

        if rid not in thresholds:
            raise KeyError(
                f"Selected product {rid} is absent from Stage-B thresholds."
            )

        maximum = float(
            products[rid]
        )

        if maximum <= tolerance:
            raise ValueError(
                f"Selected product {rid} has invalid M_i={maximum}."
            )

        q_name = (
            f"stageC_q_{index}_"
            f"{_safe_name(rid)}"
        )

        community.solver.add_variable(
            q_name,
            0,
            1,
            vartype=VarType.CONTINUOUS,
        )

        normalized_variables[rid] = q_name

    community.solver.update()

    # Fixed-z=1 forms of the Stage-A and Stage-B product constraints:
    #
    #     v_i >= alpha M_i
    #     v_i >= M_i q_i
    for index, rid in enumerate(selected_products):

        maximum = float(
            products[rid]
        )

        threshold = float(
            thresholds[rid]
        )

        q_name = normalized_variables[rid]

        community.solver.add_constraint(
            (
                f"c_stageC_product_threshold_"
                f"{index}_{_safe_name(rid)}"
            ),
            {
                rid: 1,
            },
            ">",
            threshold,
        )

        community.solver.add_constraint(
            (
                f"c_stageC_q_link_"
                f"{index}_{_safe_name(rid)}"
            ),
            {
                rid: 1,
                q_name: -maximum,
            },
            ">",
            0,
        )

    q_floor = max(
        0.0,
        q_star - lexicographic_tolerance,
    )

    community.solver.add_constraint(
        "c_stageC_preserve_stageB_objective",
        {
            q_name: 1
            for q_name in normalized_variables.values()
        },
        ">",
        q_floor,
    )

    community.solver.update()

    product_and_q_values = (
        selected_products
        + list(normalized_variables.values())
    )

    logging.info(
        "Stage C fresh model: %i selected products.",
        len(selected_products),
    )

    logging.info(
        "Stage C fresh model: sum(q_i) >= %.12g.",
        q_floor,
    )

    # --------------------------------------------------------------
    # C0: fresh biological model before split variables.
    # --------------------------------------------------------------
    if check_feasibility:
        print(
            "Stage C.1: checking feasibility before split-flux variables...",
            flush=True,
        )
        logging.info(
            "Stage C0: checking fresh biological model before split variables."
        )

        c0_solution = community.solver.solve(
            objective={},
            get_values=product_and_q_values,
            minimize=True,
        )

        if c0_solution.status != Status.OPTIMAL:
            _raise_with_iis(
                community=community,
                stage_name="Stage C0 fresh model",
                filename="stage_c0_fresh_iis.ilp",
                status=c0_solution.status,
            )

        logging.info(
            "Stage C0 fresh model is feasible."
        )
        print(
            "Stage C.1 COMPLETE: fresh biological model is feasible.",
            flush=True,
        )

        # Diagnostic repeat with no model changes. If this second solve fails,
        # the solver is mutating state simply by solving C0.
        print(
            "Stage C.1b: repeating the identical C0 feasibility solve...",
            flush=True,
        )
        c0_repeat_solution = community.solver.solve(
            objective={},
            get_values=product_and_q_values,
            minimize=True,
        )
        print(
            f"Stage C.1b status: {c0_repeat_solution.status}.",
            flush=True,
        )
        if c0_repeat_solution.status != Status.OPTIMAL:
            _raise_with_iis(
                community=community,
                stage_name="Stage C0 repeat without model changes",
                filename="stage_c0_repeat_iis.ilp",
                status=c0_repeat_solution.status,
            )

    # --------------------------------------------------------------
    # Add absolute-value epigraph variables.
    #
    # For each biological flux v_r introduce a_r >= 0 with:
    #
    #     a_r >=  v_r
    #     a_r >= -v_r
    #
    # Minimizing sum(a_r) then gives a_r = |v_r| at optimum.
    # Unlike an exact split-flux equality, these auxiliary constraints
    # cannot remove any feasible value of the original flux v_r.
    # --------------------------------------------------------------
    reaction_ids = sorted(
        community.merged_model.reactions.keys()
    )

    absolute_variables = {}

    logging.info(
        "Stage C: creating absolute-flux variables for %i reactions.",
        len(reaction_ids),
    )

    print(
        f"Stage C.2: creating 1 absolute-flux variable for each of "
        f"{len(reaction_ids)} reactions...",
        flush=True,
    )

    for index, rid in enumerate(reaction_ids):

        absolute_name = f"pfba_abs_{index}"

        community.solver.add_variable(
            absolute_name,
            0,
            float("inf"),
            vartype=VarType.CONTINUOUS,
        )

        absolute_variables[rid] = absolute_name

        completed = index + 1
        if completed % 50000 == 0 or completed == len(reaction_ids):
            print(
                f"Stage C.2 absolute variables: "
                f"{completed}/{len(reaction_ids)} reactions.",
                flush=True,
            )

    community.solver.update()

    # Diagnostic checkpoint: variables only, no epigraph constraints yet.
    if check_feasibility:
        print(
            "Stage C.2b: checking feasibility after adding auxiliary variables only...",
            flush=True,
        )
        variables_only_solution = community.solver.solve(
            objective={},
            get_values=product_and_q_values,
            minimize=True,
        )
        print(
            f"Stage C.2b status: {variables_only_solution.status}.",
            flush=True,
        )
        if variables_only_solution.status != Status.OPTIMAL:
            _raise_with_iis(
                community=community,
                stage_name="Stage C variables-only checkpoint",
                filename="stage_c_variables_only_iis.ilp",
                status=variables_only_solution.status,
            )

    print(
        "Stage C.3: adding absolute-value epigraph constraints in diagnostic batches...",
        flush=True,
    )

    # Pinpoint diagnostic:
    # - first 250 reactions are known feasible as a cumulative block;
    # - from reaction 251 onward, solve after every newly-added pair of
    #   epigraph constraints until the first infeasible reaction is found.
    for index, rid in enumerate(reaction_ids):

        absolute_name = absolute_variables[rid]

        # a_r >= v_r  ->  a_r - v_r >= 0
        community.solver.add_constraint(
            f"c_pfba_abs_pos_{index}",
            {
                absolute_name: 1,
                rid: -1,
            },
            ">",
            0,
        )

        # a_r >= -v_r  ->  a_r + v_r >= 0
        community.solver.add_constraint(
            f"c_pfba_abs_neg_{index}",
            {
                absolute_name: 1,
                rid: 1,
            },
            ">",
            0,
        )

        completed = index + 1

        # Check once at 250, then every reaction from 251 through 500.
        # If all 251-500 remain feasible, continue every 25 reactions
        # thereafter to avoid an excessive number of solves.
        should_check = (
            completed == 250
            or (251 <= completed <= 500)
            or (completed > 500 and completed % 25 == 0)
            or completed == len(reaction_ids)
        )

        if should_check:
            community.solver.update()
            print(
                f"Stage C.3 pinpoint: constraints added through "
                f"reaction {completed}/{len(reaction_ids)} "
                f"(index={index}, rid={rid}). Checking feasibility...",
                flush=True,
            )

            if check_feasibility:
                pinpoint_solution = community.solver.solve(
                    objective={},
                    get_values=False,
                    minimize=True,
                )
                print(
                    f"Stage C.3 pinpoint {completed} status: "
                    f"{pinpoint_solution.status}.",
                    flush=True,
                )

                if pinpoint_solution.status != Status.OPTIMAL:
                    previous_rid = (
                        reaction_ids[index - 1]
                        if index > 0
                        else None
                    )
                    print(
                        "FIRST FAILING REACTION IDENTIFIED:",
                        flush=True,
                    )
                    print(
                        f"  index={index}",
                        flush=True,
                    )
                    print(
                        f"  reaction_number={completed}",
                        flush=True,
                    )
                    print(
                        f"  rid={rid}",
                        flush=True,
                    )
                    print(
                        f"  previous_rid={previous_rid}",
                        flush=True,
                    )
                    _raise_with_iis(
                        community=community,
                        stage_name=(
                            "Stage C absolute-flux pinpoint at "
                            f"reaction {completed}/{len(reaction_ids)} "
                            f"({rid})"
                        ),
                        filename=(
                            f"stage_c_abs_pinpoint_{completed}_"
                            f"{_safe_name(rid)}_iis.ilp"
                        ),
                        status=pinpoint_solution.status,
                    )

    community.solver.update()

    # --------------------------------------------------------------
    # C1: epigraph representation added, still no pFBA objective.
    # --------------------------------------------------------------
    if check_feasibility:
        print(
            "Stage C.4: checking feasibility of the absolute-flux model...",
            flush=True,
        )
        logging.info(
            "Stage C1: checking feasibility after absolute-flux constraints."
        )

        c1_solution = community.solver.solve(
            objective={},
            get_values=product_and_q_values,
            minimize=True,
        )

        if c1_solution.status != Status.OPTIMAL:
            _raise_with_iis(
                community=community,
                stage_name="Stage C1 fresh absolute-flux model",
                filename="stage_c1_fresh_absolute_iis.ilp",
                status=c1_solution.status,
            )

        logging.info(
            "Stage C1 fresh absolute-flux model is feasible."
        )
        print(
            "Stage C.4 COMPLETE: absolute-flux model is feasible.",
            flush=True,
        )

    # --------------------------------------------------------------
    # C2: actual pFBA objective.
    # --------------------------------------------------------------
    pfba_objective = {
        absolute_variables[rid]: 1
        for rid in reaction_ids
    }

    logging.info(
        "Stage C2: minimizing total absolute flux."
    )

    print(
        f"Stage C.5: solving pFBA over {len(reaction_ids)} reactions...",
        flush=True,
    )

    solution = community.solver.solve(
        objective=pfba_objective,
        get_values=product_and_q_values,
        minimize=True,
    )

    if solution.status != Status.OPTIMAL:
        _raise_with_iis(
            community=community,
            stage_name="Stage C2 fresh pFBA",
            filename="stage_c2_fresh_pfba_iis.ilp",
            status=solution.status,
        )

    stage_c_q_sum = sum(
        solution.values.get(
            q_name,
            0.0,
        )
        for q_name in normalized_variables.values()
    )

    if stage_c_q_sum + tolerance < q_floor:
        raise RuntimeError(
            "Stage C solution does not preserve the Stage-B objective. "
            f"Observed sum(q_i)={stage_c_q_sum}, "
            f"required>={q_floor}."
        )

    reference_products = {}

    for rid in selected_products:

        flux = float(
            solution.values.get(
                rid,
                0.0,
            )
        )

        q_value = float(
            solution.values.get(
                normalized_variables[rid],
                0.0,
            )
        )

        reference_products[rid] = {
            "max_individual": float(
                products[rid]
            ),
            "required_fraction": float(
                stage_b["fraction"]
            ),
            "required_flux": float(
                thresholds[rid]
            ),
            "reference_flux": flux,
            "normalized_flux": (
                flux
                / float(products[rid])
            ),
            "q_value": q_value,
        }

    pfba_objective_value = float(
        solution.fobj
    )

    logging.info(
        "Stage C2 pFBA objective: %.12g",
        pfba_objective_value,
    )

    logging.info(
        "Stage C2 retained sum(q_i): %.12g",
        stage_c_q_sum,
    )

    print(
        f"Stage C COMPLETE: pFBA objective={pfba_objective_value:.12g}; "
        f"retained sum(q_i)={stage_c_q_sum:.12g}; "
        f"products={len(selected_products)}.",
        flush=True,
    )

    return {
        "selected_products": selected_products,
        "reference_products": reference_products,
        "stage_b_q_star": q_star,
        "lexicographic_tolerance": lexicographic_tolerance,
        "stage_b_objective_floor": q_floor,
        "stage_c_q_sum": stage_c_q_sum,
        "pfba_objective": pfba_objective_value,
        "exchange_fluxes": {
            rid: float(
                solution.values.get(
                    rid,
                    0.0,
                )
            )
            for rid in selected_products
        },
        "stage_c_solution": solution,
    }
