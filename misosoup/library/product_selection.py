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
    reference flux distribution using the classical split-flux pFBA
    formulation:

        v_r = v_r^+ - v_r^-

        v_r^+ >= 0
        v_r^- >= 0

        min sum_r (v_r^+ + v_r^-)

    At the optimum:

        v_r^+ + v_r^- = |v_r|

Important
---------
These functions mutate the solver contained in ``community``. They must
therefore be called sequentially, in the order A -> B.1 -> B.2 -> C, on
the same fresh LayeredCommunity object.
"""

import logging
import math
import re
from pathlib import Path

from gurobipy import GRB
from reframed.solvers.solution import Status
from reframed.solvers.solver import VarType

from ..reframed.layered_community import LayeredCommunity


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
    lexicographic_tolerance: float = 1e-5,
    tolerance: float = 1e-6,
    check_feasibility: bool = True,
) -> dict:
    """Stage C: fix the Stage-B product set and solve pFBA.

    Stage B solved the secondary optimization problem:

        Q* = max sum_i q_i

    subject to the Stage-A product constraints and:

        sum_i z_i = K*
        q_i <= z_i
        z_i = 1  =>  v_i >= M_i q_i

    Stage C fixes the identities selected in Stage B and preserves the
    Stage-B objective lexicographically:

        sum_i q_i >= Q* - epsilon

    where:

        epsilon = lexicographic_tolerance

    This is intentionally different from reconstructing the Stage-B
    objective as:

        sum_i v_i / M_i

    because Stage B optimized q_i, not v_i / M_i. Reusing the same q_i
    variables preserves exactly the objective that was optimized in Stage B.

    Importantly, Stage C does NOT add a second explicit constraint

        v_i >= alpha M_i

    for selected products. That condition is already present through the
    Stage-A indicator constraint:

        z_i = 1  =>  v_i >= alpha M_i

    and z_i is fixed to its Stage-B value in Stage C. Avoiding a duplicate
    linear version of this constraint prevents us from accidentally making
    Stage C numerically stricter than Stage B.

    Classical split-flux pFBA
    --------------------------

    For every metabolic reaction r introduce two non-negative variables:

        v_r^+ >= 0
        v_r^- >= 0

    linked to the original reaction flux by:

        v_r = v_r^+ - v_r^-

    or equivalently:

        v_r - v_r^+ + v_r^- = 0

    Stage C then solves:

        min sum_r (v_r^+ + v_r^-)

    At optimum:

        v_r^+ + v_r^- = |v_r|

    so this is equivalent to:

        min sum_r |v_r|

    Parameters
    ----------
    community
        The same LayeredCommunity used by Stages A and B. It already
        contains the Stage-A and Stage-B variables and constraints.

    stage_b
        State dictionary returned by ``maximizeProducts``.

    product_selection
        Dictionary returned by
        ``getSelectedProductsFromProductMaximization``.

    lexicographic_tolerance
        epsilon in:

            sum_i q_i >= Q* - epsilon

        Default: 1e-5.

    tolerance
        Numerical tolerance used for final validation.

    check_feasibility
        If True, solve a feasibility problem after fixing the products and
        imposing the lexicographic floor, but before adding pFBA variables.

    Returns
    -------
    dict
        Final product reference fluxes and pFBA diagnostics.
    """

    if lexicographic_tolerance < 0:
        raise ValueError(
            "lexicographic_tolerance must be >= 0."
        )

    products = stage_b["products"]
    product_ids = stage_b["product_ids"]
    product_variables = stage_b["product_variables"]
    normalized_variables = stage_b["normalized_variables"]
    q_star = stage_b["q_star"]
    k_star = stage_b["k_star"]

    selected_products = product_selection[
        "selected_products"
    ]

    selected_set = set(selected_products)

    if len(selected_set) != k_star:
        raise RuntimeError(
            "Stage C received an unexpected product count: "
            f"{len(selected_set)} instead of K*={k_star}."
        )

    # ==================================================================
    # C.1 FIX THE PRODUCT IDENTITIES SELECTED BY STAGE B
    # ==================================================================
    #
    # Stage B was allowed to choose which K* products were active.
    # Stage C must not change that biological decision.
    #
    # Therefore:
    #
    #     z_i = 1   for products selected in Stage B
    #     z_i = 0   for all other candidate products
    #
    # The original Stage-A indicator constraints remain active, so fixing
    # z_i = 1 automatically preserves:
    #
    #     v_i >= alpha M_i
    #
    # without adding a numerically different duplicate constraint.
    # ==================================================================

    for index, rid in enumerate(product_ids):
        z_name = product_variables[rid]

        fixed_value = (
            1
            if rid in selected_set
            else 0
        )

        community.solver.add_constraint(
            f"c_stageC_fix_z_{index}_{_safe_name(rid)}",
            {
                z_name: 1,
            },
            "=",
            fixed_value,
        )

    community.solver.update()

    # ==================================================================
    # C.2 LEXICOGRAPHIC RETENTION OF THE STAGE-B OBJECTIVE
    # ==================================================================
    #
    # Stage B optimized:
    #
    #     Q* = max sum_i q_i
    #
    # Therefore Stage C preserves exactly that objective, allowing only
    # a small numerical degradation epsilon:
    #
    #     sum_i q_i >= Q* - epsilon
    #
    # This is the standard lexicographic idea:
    #
    #   1. optimize the higher-priority objective,
    #   2. preserve it within a numerical tolerance,
    #   3. optimize the lower-priority objective.
    #
    # Since Stage B already contains:
    #
    #     q_i <= z_i
    #
    # and all unselected products now have z_i = 0, their q_i variables
    # are automatically forced to zero. Consequently, summing over all
    # q_i is exactly the same Stage-B objective.
    # ==================================================================

    q_floor = max(
        0.0,
        q_star - lexicographic_tolerance,
    )

    q_expression = {
        q_name: 1
        for q_name in normalized_variables.values()
    }

    community.solver.add_constraint(
        "c_stageC_preserve_stageB_objective",
        q_expression,
        ">",
        q_floor,
    )

    community.solver.update()

    logging.info(
        "Stage C: fixed %i selected products.",
        len(selected_products),
    )

    logging.info(
        "Stage C: preserving Stage-B objective with epsilon = %.12g.",
        lexicographic_tolerance,
    )

    logging.info(
        "Stage C: sum(q_i) >= Q* - epsilon = %.12g.",
        q_floor,
    )

    # ==================================================================
    # C.3 C0 FEASIBILITY CHECK
    # ==================================================================
    #
    # Before adding any pFBA variables, verify that:
    #
    #     - the Stage-B product identities are fixed, and
    #     - sum_i q_i >= Q* - epsilon
    #
    # remains feasible.
    #
    # If C0 fails, the problem is unrelated to the split-flux pFBA
    # representation.
    # ==================================================================

    if check_feasibility:
        logging.info(
            "Stage C0: checking feasibility before adding pFBA variables."
        )

        c0_solution = community.solver.solve(
            objective={},
            get_values=False,
            minimize=True,
        )

        if c0_solution.status != Status.OPTIMAL:
            _raise_with_iis(
                community=community,
                stage_name="Stage C0",
                filename="stage_c0_iis.ilp",
                status=c0_solution.status,
            )

        logging.info(
            "Stage C0 is feasible."
        )

    # ==================================================================
    # C.4 CREATE CLASSICAL SPLIT-FLUX VARIABLES
    # ==================================================================
    #
    # For every reaction r:
    #
    #     v_r^+ >= 0
    #     v_r^- >= 0
    #
    # These variables represent the forward and reverse contributions to
    # the net reaction flux.
    # ==================================================================

    reaction_ids = sorted(
        community.merged_model.reactions.keys()
    )

    positive_variables = {}
    negative_variables = {}

    logging.info(
        "Stage C: creating split-flux variables for %i reactions.",
        len(reaction_ids),
    )

    for index, rid in enumerate(reaction_ids):
        positive_name = f"pfba_pos_{index}"
        negative_name = f"pfba_neg_{index}"

        community.solver.add_variable(
            positive_name,
            0,
            math.inf,
            vartype=VarType.CONTINUOUS,
        )

        community.solver.add_variable(
            negative_name,
            0,
            math.inf,
            vartype=VarType.CONTINUOUS,
        )

        positive_variables[rid] = positive_name
        negative_variables[rid] = negative_name

    community.solver.update()

    # ==================================================================
    # C.5 LINK NET FLUX TO THE TWO NON-NEGATIVE VARIABLES
    # ==================================================================
    #
    # For every reaction:
    #
    #     v_r = v_r^+ - v_r^-
    #
    # written as the linear equality:
    #
    #     v_r - v_r^+ + v_r^- = 0
    # ==================================================================

    for index, rid in enumerate(reaction_ids):
        positive_name = positive_variables[rid]
        negative_name = negative_variables[rid]

        community.solver.add_constraint(
            f"c_pfba_split_flux_{index}",
            {
                rid: 1,
                positive_name: -1,
                negative_name: 1,
            },
            "=",
            0,
        )

    community.solver.update()

    # ==================================================================
    # C.6 pFBA OBJECTIVE
    # ==================================================================
    #
    # Solve:
    #
    #     min sum_r (v_r^+ + v_r^-)
    #
    # Because both split variables have positive objective coefficients,
    # the optimizer has no incentive to make both non-zero simultaneously.
    #
    # Therefore, at optimum:
    #
    #     v_r^+ + v_r^- = |v_r|
    #
    # and this objective is equivalent to:
    #
    #     min sum_r |v_r|
    # ==================================================================

    pfba_objective = {}

    for rid in reaction_ids:
        pfba_objective[
            positive_variables[rid]
        ] = 1

        pfba_objective[
            negative_variables[rid]
        ] = 1

    # Retrieve product fluxes, z_i values, and q_i values so that we can
    # validate the final lexicographic solution explicitly.
    values_to_get = (
        product_ids
        + list(product_variables.values())
        + list(normalized_variables.values())
    )

    logging.info(
        "Stage C: solving classical split-flux pFBA."
    )

    solution = community.solver.solve(
        objective=pfba_objective,
        get_values=values_to_get,
        minimize=True,
    )

    if solution.status != Status.OPTIMAL:
        _raise_with_iis(
            community=community,
            stage_name="Stage C pFBA",
            filename="stage_c_pfba_iis.ilp",
            status=solution.status,
        )

    # ==================================================================
    # C.7 VALIDATE THE STAGE-B OBJECTIVE AFTER pFBA
    # ==================================================================
    #
    # Reconstruct:
    #
    #     Q_C = sum_i q_i
    #
    # and verify:
    #
    #     Q_C >= Q* - epsilon
    # ==================================================================

    stage_c_q_sum = sum(
        solution.values.get(
            q_name,
            0.0,
        )
        for q_name in normalized_variables.values()
    )

    if (
        stage_c_q_sum
        + tolerance
        < q_floor
    ):
        raise RuntimeError(
            "Stage C solution does not preserve the Stage-B objective. "
            f"Observed sum(q_i)={stage_c_q_sum}, "
            f"required>={q_floor}."
        )

    # ==================================================================
    # C.8 EXTRACT FINAL REFERENCE PRODUCT FLUXES
    # ==================================================================

    reference_products = {}

    for rid in selected_products:
        flux = solution.values.get(
            rid,
            0.0,
        )

        q_value = solution.values.get(
            normalized_variables[rid],
            0.0,
        )

        reference_products[rid] = {
            "max_individual": products[rid],
            "required_fraction": stage_b["fraction"],
            "required_flux": stage_b["thresholds"][rid],
            "reference_flux": flux,
            "normalized_flux": (
                flux / products[rid]
            ),
            "q_value": q_value,
        }

    pfba_objective_value = float(
        solution.fobj
    )

    logging.info(
        "Stage C pFBA objective: %.12g",
        pfba_objective_value,
    )

    logging.info(
        "Stage C retained sum(q_i): %.12g",
        stage_c_q_sum,
    )

    exchange_fluxes = {
        rid: solution.values.get(
            rid,
            0.0,
        )
        for rid in product_ids
    }

    return {
        "selected_products": selected_products,
        "reference_products": reference_products,
        "stage_b_q_star": q_star,
        "lexicographic_tolerance": lexicographic_tolerance,
        "stage_b_objective_floor": q_floor,
        "stage_c_q_sum": stage_c_q_sum,
        "pfba_objective": pfba_objective_value,
        "exchange_fluxes": exchange_fluxes,
        "stage_c_solution": solution,
    }
