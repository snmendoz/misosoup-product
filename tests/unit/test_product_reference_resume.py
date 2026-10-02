"""Unit tests for resumable product-reference exchange scans."""

from types import SimpleNamespace

from reframed.solvers.solution import Status

from misosoup.library.product_reference import find_producible_exchanges


class _Problem:
    NumIntVars = 0
    NumBinVars = 0
    NumVars = 2
    NumConstrs = 1
    NumGenConstrs = 0

    def update(self):
        return None


class _Solver:
    def __init__(self):
        self.problem = _Problem()
        self.calls = []

    def solve(self, objective, get_values, minimize):
        rid = next(iter(objective))
        self.calls.append(rid)
        return SimpleNamespace(
            status=Status.OPTIMAL,
            values={rid: 3.0},
        )


def _community():
    solver = _Solver()
    community = SimpleNamespace(
        merged_model=SimpleNamespace(
            reactions={
                "R_EX_a_e0": object(),
                "R_EX_b_e0": object(),
            }
        ),
        solver=solver,
    )
    return community, solver


def test_product_scan_reuses_completed_exchange_checkpoint():
    """A completed exchange must not be re-solved after resume."""
    community, solver = _community()
    callback_records = {}

    result = find_producible_exchanges(
        community,
        tolerance=1e-6,
        exchange_reactions=[
            "R_EX_a_e0",
            "R_EX_b_e0",
        ],
        completed_results={
            "R_EX_a_e0": {
                "status": "optimal",
                "maximum": 2.0,
                "solve_seconds": 10.0,
            }
        },
        progress_callback=lambda rid, record: callback_records.update(
            {rid: record}
        ),
    )

    assert result == {
        "R_EX_a_e0": 2.0,
        "R_EX_b_e0": 3.0,
    }
    assert solver.calls == ["R_EX_b_e0"]
    assert callback_records["R_EX_b_e0"]["status"] == "optimal"
    assert callback_records["R_EX_b_e0"]["maximum"] == 3.0


def test_product_scan_retries_nonoptimal_checkpoint_record():
    """Non-optimal prior attempts must be retried rather than reused."""
    community, solver = _community()

    result = find_producible_exchanges(
        community,
        tolerance=1e-6,
        exchange_reactions=["R_EX_a_e0"],
        completed_results={
            "R_EX_a_e0": {
                "status": "nonoptimal",
                "maximum": None,
                "solve_seconds": 5.0,
            }
        },
    )

    assert result == {"R_EX_a_e0": 3.0}
    assert solver.calls == ["R_EX_a_e0"]
