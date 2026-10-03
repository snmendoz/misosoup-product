"""Layered microbial community implemented with COBRApy + direct Gurobi."""

import logging
import math

from cobra import Model, Reaction, Metabolite
from gurobipy import Env

from .adapters import CobraModelAdapter
from .solver import GurobiCobraSolver, VarType

BOUND_INF = 1000


def _copy_annotations(source, target):
    try:
        target.annotation = dict(source.annotation or {})
    except Exception:
        pass


class LayeredCommunity:
    """Community model with organism-local extracellular layers."""

    default_environment = None

    def __init__(
        self,
        community_id,
        models,
        env: Env = None,
        copy_models=False,
        suffix="_i",
        params=None,
        create_solver=True,
    ):
        del copy_models
        self.id = community_id
        self.suffix = suffix
        self.has_binary_variables = False

        adapted = [
            model if isinstance(model, CobraModelAdapter)
            else CobraModelAdapter(model)
            for model in models
        ]
        self.organisms = {model.id: model for model in adapted}
        if len(self.organisms) != len(adapted):
            raise ValueError("Community organism IDs must be unique.")

        self.reaction_map = {}
        self.metabolite_map = {}
        self._cobra_merged_model = self._merge_models()
        self.merged_model = CobraModelAdapter(
            self._cobra_merged_model,
            biomass_reaction="community_growth",
        )

        if create_solver:
            if env is None:
                if type(self).default_environment is None:
                    # Configure the environment before it starts so Gurobi
                    # does not emit license/banner text to stdout. MiSoSoup's
                    # CLI reserves stdout for machine-readable YAML.
                    quiet_env = Env(empty=True)
                    quiet_env.setParam("OutputFlag", 0)
                    quiet_env.setParam("LogToConsole", 0)
                    quiet_env.setParam("Method", 1)
                    quiet_env.start()
                    type(self).default_environment = quiet_env
                env = type(self).default_environment
            self.solver = GurobiCobraSolver(
                self._cobra_merged_model,
                env=env,
                params=params,
            )
        else:
            self.solver = None

    def _merge_models(self):
        community = Model(self.id)

        biomass_met = Metabolite(
            "community_biomass",
            name="Total community biomass",
            compartment="ext",
        )
        community.add_metabolites([biomass_met])

        growth = Reaction("community_growth")
        growth.name = "Community growth rate"
        growth.lower_bound = 0
        growth.upper_bound = BOUND_INF
        growth.add_metabolites({biomass_met: -1})
        community.add_reactions([growth])

        shared_external = {}

        for org_id, adapter in self.organisms.items():
            source = adapter.cobra_model
            external_compartments = {
                cid
                for cid, compartment in adapter.compartments.items()
                if compartment.external
            }

            local_metabolites = {}
            for metabolite in source.metabolites:
                local_id = f"{metabolite.id}_{org_id}"
                local = Metabolite(
                    local_id,
                    formula=metabolite.formula,
                    name=metabolite.name,
                    charge=metabolite.charge,
                    compartment=f"{metabolite.compartment}_{org_id}",
                )
                _copy_annotations(metabolite, local)
                local_metabolites[metabolite.id] = local
                self.metabolite_map[(org_id, metabolite.id)] = local_id

                if metabolite.compartment in external_compartments:
                    if metabolite.id not in shared_external:
                        shared = Metabolite(
                            metabolite.id,
                            formula=metabolite.formula,
                            name=metabolite.name,
                            charge=metabolite.charge,
                            compartment="ext",
                        )
                        _copy_annotations(metabolite, shared)
                        shared_external[metabolite.id] = shared

            community.add_metabolites(list(local_metabolites.values()))

            exchange_ids = {reaction.id for reaction in source.exchanges}

            for reaction in source.reactions:
                new_id = f"{reaction.id}_{org_id}"
                is_exchange = (
                    reaction.id.startswith("R_EX")
                    or reaction.id in exchange_ids
                )

                if is_exchange:
                    new_id += self.suffix
                    if len(reaction.metabolites) != 1:
                        raise ValueError(
                            f"Exchange {reaction.id} in {org_id} does not "
                            "have exactly one metabolite."
                        )
                    original_met = next(iter(reaction.metabolites))
                    shared = shared_external[original_met.id]
                    local = local_metabolites[original_met.id]

                    new_reaction = Reaction(new_id)
                    new_reaction.lower_bound = -BOUND_INF
                    new_reaction.upper_bound = BOUND_INF
                    new_reaction.add_metabolites({shared: 1, local: -1})
                else:
                    new_reaction = Reaction(new_id)
                    new_reaction.name = reaction.name
                    new_reaction.lower_bound = float(reaction.lower_bound)
                    new_reaction.upper_bound = float(reaction.upper_bound)
                    new_reaction.add_metabolites({
                        local_metabolites[met.id]: float(coefficient)
                        for met, coefficient in reaction.metabolites.items()
                    })

                    if reaction.id == adapter.biomass_reaction:
                        new_reaction.add_metabolites({biomass_met: 1})

                    try:
                        new_reaction.gene_reaction_rule = (
                            reaction.gene_reaction_rule
                        )
                    except Exception:
                        pass

                _copy_annotations(reaction, new_reaction)
                community.add_reactions([new_reaction])
                self.reaction_map[(org_id, reaction.id)] = new_id

        existing_metabolites = {met.id for met in community.metabolites}
        community.add_metabolites([
            met
            for mid, met in shared_external.items()
            if mid not in existing_metabolites
        ])

        existing_reactions = {reaction.id for reaction in community.reactions}
        for metabolite_id, metabolite in shared_external.items():
            rid = (
                f"R_EX_{metabolite_id[2:]}"
                if metabolite_id.startswith("M_")
                else f"R_EX_{metabolite_id}"
            )
            if rid in existing_reactions:
                continue
            exchange = Reaction(rid)
            exchange.lower_bound = -BOUND_INF
            exchange.upper_bound = BOUND_INF
            exchange.add_metabolites({metabolite: -1})
            community.add_reactions([exchange])

        return community

    def setup_binary_variables(self, minimal_growth):
        for org_id in self.organisms:
            self.solver.add_variable(
                f"y_{org_id}",
                0,
                1,
                vartype=VarType.BINARY,
            )
        self.solver.update()

        for org_id, org_model in self.organisms.items():
            org_var = f"y_{org_id}"
            for r_id, reaction in org_model.reactions.items():
                if (
                    not r_id.startswith("R_EX")
                    and r_id != org_model.biomass_reaction
                    and reaction.lb * reaction.ub <= 0
                ):
                    continue

                merged_id = self.reaction_map[(org_id, r_id)]
                upper = BOUND_INF
                lower = -BOUND_INF

                if r_id == org_model.biomass_reaction:
                    lower = minimal_growth

                if reaction.lb * reaction.ub > 0:
                    lower = (
                        -BOUND_INF if math.isinf(reaction.lb)
                        else reaction.lb
                    )
                    upper = (
                        BOUND_INF if math.isinf(reaction.ub)
                        else reaction.ub
                    )
                    self.solver.set_bounds({
                        merged_id: (-BOUND_INF, BOUND_INF)
                    })

                self.solver.add_constraint(
                    f"c_{merged_id}_lb",
                    {merged_id: 1, org_var: -lower},
                    ">",
                    0,
                )
                self.solver.add_constraint(
                    f"c_{merged_id}_ub",
                    {merged_id: 1, org_var: -upper},
                    "<",
                    0,
                )

        self.solver.update()
        self.has_binary_variables = True

    def setup_growth_requirement(self, minimal_growth):
        for org_id, org_model in self.organisms.items():
            merged_id = self.reaction_map[
                (org_id, org_model.biomass_reaction)
            ]
            self.solver.add_constraint(
                f"c_{merged_id}_lb",
                {merged_id: 1},
                ">",
                minimal_growth,
            )
        self.solver.update()

    def setup_medium(self, medium):
        missing = set(medium) - set(self.merged_model.reactions)
        for rid in sorted(missing):
            logging.warning("Missing reaction %s in model %s", rid, self.id)

        for rid in self.merged_model.reactions:
            if rid.startswith("R_EX_") and not rid.endswith("_i"):
                bound = medium[rid] if rid in medium else 0
                self.solver.add_constraint(
                    f"c_{rid}_lb",
                    {rid: 1},
                    ">",
                    bound,
                )
        self.solver.update()

    def setup_parsimony(self):
        for rid in self.merged_model.reactions:
            self.solver.add_variable(f"abs_{rid}", 0, math.inf)
        self.solver.update()

        for rid in self.merged_model.reactions:
            self.solver.add_constraint(
                f"c_{rid}_abs_pos",
                {f"abs_{rid}": 1, rid: -1},
                ">",
                0,
            )
            self.solver.add_constraint(
                f"c_{rid}_abs_neg",
                {f"abs_{rid}": 1, rid: 1},
                ">",
                0,
            )
        self.solver.update()

    def check_feasibility(self, values):
        existing = set(values) & set(self.merged_model.reactions)
        return self.solver.solve(get_values=existing)

    def objective_optimization(self, objective, values):
        existing = set(values) & set(self.merged_model.reactions)
        return self.solver.solve(
            objective=objective,
            get_values=existing,
            minimize=False,
        )

    def parsimony_optimization(self, objective, value, values):
        if value > 0:
            self.solver.add_constraint("c_growth", objective, ">", value)
            self.solver.update()

        existing = set(values) & set(self.merged_model.reactions)
        solution = self.solver.solve(
            objective={
                f"abs_{rid}": 1
                for rid in self.merged_model.reactions
            },
            get_values=existing,
            minimize=True,
        )
        self.solver.remove_constraint("c_growth")
        return solution
