"""Unit tests for biological prefiltering of product exchanges."""

from types import SimpleNamespace

from misosoup.library.product_filter import filter_exchange_candidates


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
