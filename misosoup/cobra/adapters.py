"""Adapters exposing the historical MiSoSoup model API on COBRApy."""

from dataclasses import dataclass

from cobra.util.solver import linear_reaction_coefficients


@dataclass
class CompartmentAdapter:
    id: str
    name: str
    external: bool = False


class MetaboliteAdapter:
    def __init__(self, metabolite):
        self._metabolite = metabolite
        self.id = metabolite.id
        self.name = metabolite.name or metabolite.id
        self.compartment = metabolite.compartment
        self.metadata = dict(metabolite.annotation or {})
        if metabolite.formula:
            self.metadata["FORMULA"] = metabolite.formula


class ReactionAdapter:
    def __init__(self, reaction):
        self._reaction = reaction
        self.id = reaction.id
        self.name = reaction.name or reaction.id
        self.lb = float(reaction.lower_bound)
        self.ub = float(reaction.upper_bound)
        self.reversible = self.lb < 0 < self.ub
        self.stoichiometry = {
            metabolite.id: float(coefficient)
            for metabolite, coefficient in reaction.metabolites.items()
        }
        self.metadata = dict(reaction.annotation or {})
        self.gene_reaction_rule = reaction.gene_reaction_rule
        self.is_exchange = bool(
            reaction.id.startswith("R_EX") or reaction.boundary
        )


class CobraModelAdapter:
    """Mapping-style facade around a COBRApy Model."""

    def __init__(self, model, biomass_reaction=None):
        self.cobra_model = model
        self.id = model.id

        exchange_compartments = {
            metabolite.compartment
            for reaction in model.exchanges
            for metabolite in reaction.metabolites
        }

        self.compartments = {
            compartment_id: CompartmentAdapter(
                compartment_id,
                compartment_name,
                compartment_id in exchange_compartments,
            )
            for compartment_id, compartment_name in model.compartments.items()
        }
        self.metabolites = {
            metabolite.id: MetaboliteAdapter(metabolite)
            for metabolite in model.metabolites
        }
        self.reactions = {
            reaction.id: ReactionAdapter(reaction)
            for reaction in model.reactions
        }
        self.genes = {gene.id: gene for gene in model.genes}
        self.biomass_reaction = biomass_reaction or self._detect_biomass()

    def _detect_biomass(self):
        coefficients = linear_reaction_coefficients(self.cobra_model)
        nonzero = [
            reaction.id
            for reaction, coefficient in coefficients.items()
            if abs(float(coefficient)) > 0
        ]
        if len(nonzero) == 1:
            return nonzero[0]

        candidates = [
            reaction.id
            for reaction in self.cobra_model.reactions
            if "biomass" in reaction.id.lower()
            or reaction.id.lower().startswith("growth")
            or "biomass" in (reaction.name or "").lower()
        ]
        if candidates:
            return candidates[0]

        raise ValueError(
            f"Unable to identify biomass reaction for model {self.cobra_model.id}."
        )
