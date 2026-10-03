"""Regression tests for the COBRApy + direct Gurobi backend."""

import math

import pytest
from cobra import Model, Reaction, Metabolite

from misosoup.cobra.solver import GurobiCobraSolver, Status


def _toy_model():
    model = Model("toy")

    metabolite = Metabolite("M_a_c", compartment="c")

    source = Reaction("R_source")
    source.lower_bound = 0
    source.upper_bound = 10
    source.add_metabolites({metabolite: 1})

    drain = Reaction("R_drain")
    drain.lower_bound = 0
    drain.upper_bound = 1000
    drain.add_metabolites({metabolite: -1})

    model.add_reactions([source, drain])
    model.objective = drain
    return model


def test_direct_gurobi_mass_balance_and_objective():
    solver = GurobiCobraSolver(_toy_model())

    solution = solver.solve(
        objective={"R_drain": 1},
        get_values=["R_source", "R_drain"],
        minimize=False,
    )

    assert solution.status == Status.OPTIMAL
    assert solution.values["R_source"] == pytest.approx(10)
    assert solution.values["R_drain"] == pytest.approx(10)


def test_absolute_value_epigraph_preserves_feasibility():
    solver = GurobiCobraSolver(_toy_model())

    solver.add_constraint(
        "c_required_drain",
        {"R_drain": 1},
        ">",
        5,
    )
    solver.add_variable("abs_drain", 0, math.inf)
    solver.update()

    solver.add_constraint(
        "c_abs_pos",
        {"abs_drain": 1, "R_drain": -1},
        ">",
        0,
    )
    solver.add_constraint(
        "c_abs_neg",
        {"abs_drain": 1, "R_drain": 1},
        ">",
        0,
    )
    solver.update()

    feasibility = solver.solve(
        objective={},
        get_values=["R_drain", "abs_drain"],
        minimize=True,
    )
    assert feasibility.status == Status.OPTIMAL

    pfba = solver.solve(
        objective={"abs_drain": 1},
        get_values=["R_source", "R_drain", "abs_drain"],
        minimize=True,
    )

    assert pfba.status == Status.OPTIMAL
    assert pfba.values["R_drain"] == pytest.approx(5)
    assert pfba.values["R_source"] == pytest.approx(5)
    assert pfba.values["abs_drain"] == pytest.approx(5)
