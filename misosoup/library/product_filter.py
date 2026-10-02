"""Systematic filtering of candidate community products for C/N/S/P metabolism.

Biological objective
====================

The filter is designed for microbial community studies focused on the
biogeochemical/metabolic cycling of:

    C = carbon
    N = nitrogen
    S = sulfur
    P = phosphorus

A candidate metabolite is therefore retained when its molecular formula
contains at least one of these elements.

This automatically keeps ecologically relevant inorganic compounds such as:

    CO2, HCO3-, NH4+, NO2-, NO3-, H2S, sulfate, sulfite,
    elemental sulfur, phosphate, and organic C/N/S/P metabolites.

At the same time, simple nuisance species such as:

    H+, H2O, H2, Na+, K+, Fe2+, Fe3+, Cu2+, Mg2+, Ca2+, ...

are excluded because they contain none of C/N/S/P.

Molecular oxygen (O2) is treated separately because it contains none of
C/N/S/P but can be metabolically/ecologically important as an electron
acceptor or product. Its inclusion is controlled by ``keep_oxygen``.

Missing formulas are kept conservatively and explicitly flagged for review.
"""

import re
from typing import Dict, Tuple

from ..reframed.layered_community import LayeredCommunity


# Elements that define the biological scope of the product analysis.
TARGET_ELEMENTS = {"C", "N", "S", "P"}

# Molecular oxygen is deliberately treated as a separate configurable class.
OXYGEN_IDS = {
    "o2",
}


def _normalize_metabolite_id(metabolite_id: str) -> str:
    """Normalize common BiGG/ReFramed metabolite IDs.

    Examples
    --------
    M_co2_e -> co2
    co2_e   -> co2
    M_o2_e  -> o2
    """

    value = metabolite_id

    if value.startswith("M_"):
        value = value[2:]

    # BiGG/ReFramed models commonly use "_e"; gapseq models often use
    # "_e0" for the extracellular compartment. Strip the longer suffix first.
    for suffix in ("_e0", "_e"):
        if value.endswith(suffix):
            value = value[:-len(suffix)]
            break

    return value.lower()


def _parse_formula(formula: str) -> Dict[str, float]:
    """Parse a compact biochemical formula into elemental counts.

    Examples
    --------
    CO2      -> {"C": 1, "O": 2}
    H2S      -> {"H": 2, "S": 1}
    HPO4     -> {"H": 1, "P": 1, "O": 4}
    Fe       -> {"Fe": 1}

    Decimal stoichiometries are accepted for robustness.
    """

    if not formula:
        return {}

    counts = {}

    for element, raw_count in re.findall(
        r"([A-Z][a-z]?)([0-9]*\.?[0-9]*)",
        formula,
    ):
        count = 1.0 if raw_count == "" else float(raw_count)
        counts[element] = counts.get(element, 0.0) + count

    return counts


def _get_exchange_metabolite(
    community: LayeredCommunity,
    reaction_id: str,
):
    """Return the metabolite associated with a global exchange reaction.

    LayeredCommunity global exchanges have the form:

        reaction stoichiometry = {metabolite_id: -1}

    and should therefore contain exactly one metabolite.
    """

    reaction = community.merged_model.reactions[reaction_id]

    metabolite_ids = list(
        reaction.stoichiometry.keys()
    )

    if len(metabolite_ids) != 1:
        raise ValueError(
            f"Expected one metabolite in exchange {reaction_id}, "
            f"found {len(metabolite_ids)}."
        )

    metabolite_id = metabolite_ids[0]

    return community.merged_model.metabolites[
        metabolite_id
    ]


def _classify_exchange_candidate(
    community: LayeredCommunity,
    reaction_id: str,
    keep_oxygen: bool = True,
) -> dict:
    """Classify one exchange before any optimization is performed.

    This is intentionally independent of flux values so biologically irrelevant
    exchanges (for example water, protons, H2 and metals) can be removed from
    the expensive product-maximization scan.
    """
    metabolite = _get_exchange_metabolite(
        community,
        reaction_id,
    )

    metabolite_id = metabolite.id
    normalized_id = _normalize_metabolite_id(
        metabolite_id
    )
    formula = metabolite.metadata.get("FORMULA")
    element_counts = _parse_formula(formula)
    elements = set(element_counts.keys())
    matched_target_elements = sorted(
        elements & TARGET_ELEMENTS
    )

    keep = False
    reason = "outside_CNSP_scope"

    if not formula:
        keep = True
        reason = "missing_formula_kept_for_review"
    elif matched_target_elements:
        keep = True
        reason = (
            "contains_target_element:"
            + ",".join(matched_target_elements)
        )
    elif normalized_id in OXYGEN_IDS:
        keep = bool(keep_oxygen)
        reason = (
            "oxygen_kept"
            if keep_oxygen
            else "oxygen_excluded"
        )

    return {
        "metabolite_id": metabolite_id,
        "metabolite_name": metabolite.name,
        "formula": formula,
        "elements": sorted(elements),
        "target_elements": matched_target_elements,
        "keep": keep,
        "reason": reason,
    }


def filter_exchange_candidates(
    community: LayeredCommunity,
    exchange_reactions: list,
    keep_oxygen: bool = True,
) -> Tuple[list, dict]:
    """Filter global exchanges *before* individual product maximization.

    Only exchanges relevant to C/N/S/P metabolism (plus O2 when requested)
    are returned for optimization. Exchanges with missing formula are retained
    conservatively and flagged for review.
    """
    selected = []
    audit = {}

    for reaction_id in sorted(exchange_reactions):
        classification = _classify_exchange_candidate(
            community,
            reaction_id,
            keep_oxygen=keep_oxygen,
        )
        audit[reaction_id] = classification

        if classification["keep"]:
            selected.append(reaction_id)

    return selected, audit


def filter_product_candidates(
    community: LayeredCommunity,
    max_secretion: dict,
    keep_oxygen: bool = True,
) -> Tuple[dict, dict]:
    """Filter products according to C/N/S/P biological relevance.

    Mathematical definition
    -----------------------

    Let E_i be the set of elements present in metabolite i and let:

        T = {C, N, S, P}

    The default retained set is:

        L_CNSP = {
            i in L :
            E_i intersect T != empty
            OR i = O2 when keep_oxygen = True
        }

    Examples
    --------

    Retained because they contain a target element:

        CO2       -> C
        H2S       -> S
        S         -> S
        NH4       -> N
        NO3       -> N
        SO4       -> S
        PO4       -> P
        acetate   -> C
        amino acid -> C/N

    Excluded because they contain none of C/N/S/P:

        H+        -> H
        H2O       -> H/O
        H2        -> H
        Fe2+      -> Fe
        Cu2+      -> Cu
        Na+       -> Na

    O2 is configurable because it contains only oxygen but can be important
    for microbial redox metabolism.

    Missing formula policy
    ----------------------

    If the molecular formula is unavailable, the metabolite is kept and
    flagged:

        reason = "missing_formula_kept_for_review"

    This avoids silently losing biologically relevant products because of
    incomplete model annotations.
    """

    filtered_products = {}
    audit = {}

    for reaction_id, maximum in max_secretion.items():
        classification = _classify_exchange_candidate(
            community,
            reaction_id,
            keep_oxygen=keep_oxygen,
        )

        if classification["keep"]:
            filtered_products[
                reaction_id
            ] = maximum

        audit[reaction_id] = {
            **classification,
            "maximum_secretion": float(maximum),
        }

    return filtered_products, audit
