"""COBRApy/Gurobi backend with the solver API used by MiSoSoup."""

from dataclasses import dataclass
from enum import Enum
from math import inf

import gurobipy as gp
from gurobipy import GRB


class Status(Enum):
    OPTIMAL = "optimal"
    INFEASIBLE = "infeasible"
    UNBOUNDED = "unbounded"
    INF_OR_UNB = "inf_or_unbounded"
    UNKNOWN = "unknown"


class VarType(Enum):
    CONTINUOUS = "continuous"
    BINARY = "binary"
    INTEGER = "integer"


class Parameter(Enum):
    TIME_LIMIT = "time_limit"
    FEASIBILITY_TOL = "feasibility_tol"
    INT_FEASIBILITY_TOL = "int_feasibility_tol"
    OPTIMALITY_TOL = "optimality_tol"
    MIP_ABS_GAP = "mip_abs_gap"
    MIP_REL_GAP = "mip_rel_gap"


@dataclass
class Solution:
    status: Status
    fobj: float | None = None
    values: dict | None = None


_PARAM_MAP = {
    Parameter.TIME_LIMIT: GRB.Param.TimeLimit,
    Parameter.FEASIBILITY_TOL: GRB.Param.FeasibilityTol,
    Parameter.INT_FEASIBILITY_TOL: GRB.Param.IntFeasTol,
    Parameter.OPTIMALITY_TOL: GRB.Param.OptimalityTol,
    Parameter.MIP_ABS_GAP: GRB.Param.MIPGapAbs,
    Parameter.MIP_REL_GAP: GRB.Param.MIPGap,
}

_VTYPE_MAP = {
    VarType.CONTINUOUS: GRB.CONTINUOUS,
    VarType.BINARY: GRB.BINARY,
    VarType.INTEGER: GRB.INTEGER,
}


def _bound(value):
    if value == inf:
        return GRB.INFINITY
    if value == -inf:
        return -GRB.INFINITY
    return float(value)


class GurobiCobraSolver:
    """Direct Gurobi LP/MILP generated from a COBRApy model."""

    def __init__(self, cobra_model, env=None, params=None):
        self.cobra_model = cobra_model
        self.params = dict(params or {})
        self.problem = gp.Model(env=env) if env is not None else gp.Model()
        self.problem.Params.OutputFlag = 0
        self.variables = []
        self.constraints = []

        for parameter, value in self.params.items():
            grb_parameter = _PARAM_MAP.get(parameter)
            if grb_parameter is not None:
                self.problem.setParam(grb_parameter, value)

        for reaction in cobra_model.reactions:
            self.problem.addVar(
                name=reaction.id,
                lb=_bound(reaction.lower_bound),
                ub=_bound(reaction.upper_bound),
                vtype=GRB.CONTINUOUS,
            )
            self.variables.append(reaction.id)
        self.problem.update()

        for metabolite in cobra_model.metabolites:
            expression = gp.LinExpr()
            for reaction in metabolite.reactions:
                coefficient = reaction.metabolites[metabolite]
                if coefficient:
                    expression.add(
                        self.problem.getVarByName(reaction.id),
                        float(coefficient),
                    )
            self.problem.addConstr(expression == 0, name=metabolite.id)
            self.constraints.append(metabolite.id)

        self.problem.update()

    def update(self):
        self.problem.update()

    def add_variable(self, var_id, lb=0, ub=inf, vartype=VarType.CONTINUOUS):
        self.problem.update()
        if self.problem.getVarByName(var_id) is not None:
            raise ValueError(f"Variable already exists: {var_id}")
        self.problem.addVar(
            name=var_id,
            lb=_bound(lb),
            ub=_bound(ub),
            vtype=_VTYPE_MAP[vartype],
        )
        self.variables.append(var_id)

    def add_constraint(self, constr_id, lhs, sense, rhs):
        self.problem.update()
        expression = gp.LinExpr()
        for variable_id, coefficient in lhs.items():
            if not coefficient:
                continue
            variable = self.problem.getVarByName(variable_id)
            if variable is None:
                raise KeyError(
                    f"Constraint {constr_id} references unknown variable "
                    f"{variable_id}."
                )
            expression.add(variable, float(coefficient))

        if sense == "=":
            constr = self.problem.addConstr(expression == rhs, name=constr_id)
        elif sense == ">":
            constr = self.problem.addConstr(expression >= rhs, name=constr_id)
        elif sense == "<":
            constr = self.problem.addConstr(expression <= rhs, name=constr_id)
        else:
            raise ValueError(f"Unknown constraint sense: {sense}")

        self.constraints.append(constr_id)
        return constr

    def remove_constraint(self, constr_id):
        if not isinstance(constr_id, str):
            return
        self.problem.update()
        constraint = self.problem.getConstrByName(constr_id)
        if constraint is not None:
            self.problem.remove(constraint)
            if constr_id in self.constraints:
                self.constraints.remove(constr_id)
            self.problem.update()

    def set_bounds(self, bounds):
        self.problem.update()
        for variable_id, (lower, upper) in bounds.items():
            variable = self.problem.getVarByName(variable_id)
            if variable is None:
                raise KeyError(f"Unknown variable: {variable_id}")
            variable.LB = _bound(lower)
            variable.UB = _bound(upper)
        self.problem.update()

    def set_parameter(self, parameter, value):
        grb_parameter = _PARAM_MAP.get(parameter)
        if grb_parameter is None:
            raise ValueError(f"Unsupported solver parameter: {parameter}")
        self.problem.setParam(grb_parameter, value)
        self.params[parameter] = value

    def solve(self, objective=None, get_values=True, minimize=True):
        self.problem.update()
        objective = objective or {}
        if isinstance(objective, str):
            objective = {objective: 1.0}

        expression = gp.LinExpr()
        for variable_id, coefficient in objective.items():
            variable = self.problem.getVarByName(variable_id)
            if variable is None:
                raise KeyError(f"Objective references unknown variable {variable_id}.")
            if coefficient:
                expression.add(variable, float(coefficient))

        self.problem.setObjective(
            expression,
            GRB.MINIMIZE if minimize else GRB.MAXIMIZE,
        )
        self.problem.optimize()

        status = {
            GRB.OPTIMAL: Status.OPTIMAL,
            GRB.INFEASIBLE: Status.INFEASIBLE,
            GRB.UNBOUNDED: Status.UNBOUNDED,
            GRB.INF_OR_UNBD: Status.INF_OR_UNB,
        }.get(self.problem.Status, Status.UNKNOWN)

        values = None
        fobj = None
        if status == Status.OPTIMAL:
            fobj = float(self.problem.ObjVal)
            if get_values:
                names = (
                    list(get_values)
                    if isinstance(get_values, (list, tuple, set))
                    else list(self.variables)
                )
                values = {}
                for name in names:
                    variable = self.problem.getVarByName(name)
                    if variable is not None:
                        values[name] = float(variable.X)

        return Solution(status=status, fobj=fobj, values=values)
