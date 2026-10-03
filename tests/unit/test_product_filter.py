"""Unit tests for biological prefiltering of product exchanges."""

from types import SimpleNamespace

from misosoup.library.product_filter import filter_exchange_candidates


def _community(exchange_specs):
    """Build a minimal community fixture from exchange specifications."""
    metabolites = {}
    reactions = {}

    for reaction_id, (metabolite_id, name, formula) in exchange_specs.items():
        metabolites[metabolite_id] = SimpleNamespace(
            id=metabolite_id,
            name=name,
            metadata={"FORMULA": formula} if formula else {},
        )
        reactions[reaction_id] = SimpleNamespace(
            stoichiometry={metabolite_id: -1}
        )

    return SimpleNamespace(
        merged_model=SimpleNamespace(
            reactions=reactions,
            metabolites=metabolites,
        )
    )


def _community_with_exchanges():
    metabolites = {
        "M_h2o_e0": SimpleNamespace(
            id="M_h2o_e0",
            name="Water",
            metadata={"FORMULA": "H2O"},
        ),
        "M_h_e0": SimpleNamespace(
            id="M_h_e0",
            name="Proton",
            metadata={"FORMULA": "H"},
        ),
        "M_o2_e0": SimpleNamespace(
            id="M_o2_e0",
            name="Oxygen",
            metadata={"FORMULA": "O2"},
        ),
        "M_glc_e0": SimpleNamespace(
            id="M_glc_e0",
            name="Glucose",
            metadata={"FORMULA": "C6H12O6"},
        ),
        "M_unknown_e0": SimpleNamespace(
            id="M_unknown_e0",
            name="Unknown",
            metadata={},
        ),
    }

    reactions = {
        "R_EX_h2o_e0": SimpleNamespace(
            stoichiometry={"M_h2o_e0": -1}
        ),
        "R_EX_h_e0": SimpleNamespace(
            stoichiometry={"M_h_e0": -1}
        ),
        "R_EX_o2_e0": SimpleNamespace(
            stoichiometry={"M_o2_e0": -1}
        ),
        "R_EX_glc_e0": SimpleNamespace(
            stoichiometry={"M_glc_e0": -1}
        ),
        "R_EX_unknown_e0": SimpleNamespace(
            stoichiometry={"M_unknown_e0": -1}
        ),
    }

    return SimpleNamespace(
        merged_model=SimpleNamespace(
            reactions=reactions,
            metabolites=metabolites,
        )
    )


def test_prefilter_excludes_water_and_protons_before_optimization():
    """Water and protons should never reach the expensive LP product scan."""
    community = _community_with_exchanges()

    selected, audit = filter_exchange_candidates(
        community,
        list(community.merged_model.reactions),
        keep_oxygen=True,
    )

    assert "R_EX_h2o_e0" not in selected
    assert "R_EX_h_e0" not in selected
    assert audit["R_EX_h2o_e0"]["reason"] == "outside_CNSP_scope"
    assert audit["R_EX_h_e0"]["reason"] == "outside_CNSP_scope"


def test_prefilter_keeps_cnsp_oxygen_and_missing_formula():
    """C/N/S/P products, optional oxygen and unknown formulas are retained."""
    community = _community_with_exchanges()

    selected, audit = filter_exchange_candidates(
        community,
        list(community.merged_model.reactions),
        keep_oxygen=True,
    )

    assert "R_EX_glc_e0" in selected
    assert "R_EX_o2_e0" in selected
    assert "R_EX_unknown_e0" in selected
    assert audit["R_EX_o2_e0"]["reason"] == "oxygen_kept"
    assert (
        audit["R_EX_unknown_e0"]["reason"]
        == "missing_formula_kept_for_review"
    )


def test_prefilter_can_exclude_oxygen():
    """Oxygen follows the keep_oxygen switch."""
    community = _community_with_exchanges()

    selected, audit = filter_exchange_candidates(
        community,
        list(community.merged_model.reactions),
        keep_oxygen=False,
    )

    assert "R_EX_o2_e0" not in selected
    assert audit["R_EX_o2_e0"]["reason"] == "oxygen_excluded"


def test_gapseq_modelseed_oxygen_is_configurable():
    """ModelSEED cpd00007 must be treated as molecular oxygen."""
    community = _community(
        {
            "R_EX_cpd00007_e0": (
                "M_cpd00007_e0",
                "O2",
                "O2",
            ),
        }
    )

    selected, audit = filter_exchange_candidates(
        community,
        ["R_EX_cpd00007_e0"],
        keep_oxygen=True,
    )
    assert selected == ["R_EX_cpd00007_e0"]
    assert audit["R_EX_cpd00007_e0"]["reason"] == "oxygen_kept"

    selected, audit = filter_exchange_candidates(
        community,
        ["R_EX_cpd00007_e0"],
        keep_oxygen=False,
    )
    assert selected == []
    assert audit["R_EX_cpd00007_e0"]["reason"] == "oxygen_excluded"
