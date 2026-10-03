"""Regression tests for the parallel HiGHS Product Scan."""

from inspect import signature

import numpy as np
from scipy.sparse import csc_matrix

from misosoup.library.highs_product_scan import (
    ProductScanLP,
    run_parallel_product_scan,
)
from misosoup.cobra.layered_community import LayeredCommunity


def _tiny_lp():
    # x + y = 1, 0 <= x,y <= 1.
    # Therefore max(x)=1 and max(y)=1.
    reaction_ids = ["R_EX_x", "R_EX_y"]
    return ProductScanLP(
        reaction_ids=reaction_ids,
        reaction_index={rid: i for i, rid in enumerate(reaction_ids)},
        a_eq=csc_matrix([[1.0, 1.0]]),
        b_eq=np.array([1.0]),
        bounds=[(0.0, 1.0), (0.0, 1.0)],
    )


def test_layered_community_exposes_create_solver_keyword():
    params = signature(LayeredCommunity.__init__).parameters
    assert "create_solver" in params
    assert params["create_solver"].default is True
    assert hasattr(LayeredCommunity, "default_environment")


def test_highs_product_scan_single_worker():
    producible, results, _ = run_parallel_product_scan(
        _tiny_lp(),
        ["R_EX_x", "R_EX_y"],
        workers=1,
        tolerance=1e-9,
    )
    assert producible["R_EX_x"] == 1.0
    assert producible["R_EX_y"] == 1.0
    assert results["R_EX_x"]["status"] == "optimal"
    assert results["R_EX_y"]["status"] == "optimal"


def test_highs_product_scan_two_workers():
    producible, results, _ = run_parallel_product_scan(
        _tiny_lp(),
        ["R_EX_x", "R_EX_y"],
        workers=2,
        tolerance=1e-9,
    )
    assert producible["R_EX_x"] == 1.0
    assert producible["R_EX_y"] == 1.0
    assert results["R_EX_x"]["status"] == "optimal"
    assert results["R_EX_y"]["status"] == "optimal"


def test_highs_product_scan_reuses_checkpoint():
    producible, results, _ = run_parallel_product_scan(
        _tiny_lp(),
        ["R_EX_x", "R_EX_y"],
        workers=1,
        tolerance=1e-9,
        completed_results={
            "R_EX_x": {
                "status": "optimal",
                "maximum": 1.0,
                "solve_seconds": 0.1,
            }
        },
    )
    assert producible["R_EX_x"] == 1.0
    assert producible["R_EX_y"] == 1.0
    assert results["R_EX_x"]["maximum"] == 1.0
    assert results["R_EX_y"]["status"] == "optimal"
