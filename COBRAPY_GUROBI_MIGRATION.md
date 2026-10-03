# COBRApy + Gurobi backend migration

This experimental branch removes ReFramed from the active MiSoSoup runtime.

## Architecture

    SBML
      |
      v
    COBRApy
      |
      +-- model/reaction/metabolite representation
      +-- biomass/objective discovery
      +-- SBML FBC parsing
      |
      v
    MiSoSoup LayeredCommunity
      |
      v
    Gurobi direct API
      |
      +-- LP product scan
      +-- Stage A MILP
      +-- Stage B MILP
      +-- Stage C pFBA
      +-- Stage D minimum community MILP

The public LayeredCommunity and solver methods used by the existing pipeline
are intentionally preserved so that the biological formulations of Stages
A-D can be compared against the ReFramed implementation.

## Important

The old misosoup/reframed source directory is retained temporarily as a
read-only reference during validation, but it is no longer listed as an
installed package and ReFramed is no longer an install dependency.

## Validation on Leftraru

Create a clean environment or install the branch in the existing environment:

    cd ~/misosoup_product
    git fetch origin
    git switch experiment/cobrapy-gurobi-backend
    git pull --ff-only origin experiment/cobrapy-gurobi-backend

    conda activate misosoup
    python -m pip install -e .

Check the backend:

    python - <<'PY'
    import cobra
    import gurobipy
    from misosoup.cobra.layered_community import LayeredCommunity
    print("cobra", cobra.__version__)
    print("gurobi", gurobipy.gurobi.version())
    print("backend", LayeredCommunity.__module__)
    PY

Then run:

    python -m pytest -q tests
    python -m pytest -q test_product_reference.py
    python -m pytest -q test_simultaneous_products.py

## Acceptance criteria

1. No active runtime import from reframed.
2. The three marine SBML models load with the same organism IDs.
3. The full-community LP is feasible at the same minimum growth.
4. Individual product maxima agree within numerical tolerance.
5. Stage A returns the same K*.
6. Stage B returns the same selected product identities and Q* within tolerance.
7. Stage C remains feasible after pFBA variables/constraints are added.
8. Stage D returns product-preserving minimum communities.
9. The serializable output schema remains unchanged.

## Why this branch exists

The ReFramed Stage-C diagnostic showed a biologically feasible C0 solution
where R_AACPS3_I3M07 carried positive flux. Its absolute-value epigraph
constraints were algebraically correct and unbounded above, yet adding those
constraints changed Gurobi's reported status to infeasible. Rebuilding the
backend independently is therefore a useful implementation-level control.
