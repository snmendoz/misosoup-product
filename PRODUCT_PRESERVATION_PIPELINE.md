# Product-preserving minimal communities — development checkpoint

This branch extends MiSoSoup with a workflow to define a metabolic reference
phenotype for a full microbial community and then search for the
minimum-cardinality subcommunities that preserve that phenotype.

## Current workflow

1. **Reference community** — force all organisms active and require minimum
   growth for every member.
2. **Individual production maxima** — maximize each global community exchange
   independently to obtain `M_i`.
3. **Systematic C/N/S/P filter** — retain products whose molecular formula
   contains C, N, S, or P. Molecular oxygen is configurable; metabolites with
   missing formulas are retained and flagged for review.
4. **Stage A** — maximize the number of simultaneously produced products:

   `z_i = 1 => v_i >= alpha * M_i`.

5. **Stage B** — with the Stage-A cardinality fixed, maximize normalized product
   production using `q_i` variables:

   `z_i = 1 => v_i >= M_i * q_i`, then maximize `sum(q_i)`.

6. **Stage C** — rebuild the selected-product reference problem in a **fresh
   community solver**, preserve the Stage-B objective lexicographically,

   `sum(q_i) >= Q* - epsilon`,

   and run split-flux pFBA:

   `v_r = v_r+ - v_r-`, minimizing `sum(v_r+ + v_r-)`.

   Diagnostic feasibility checks are performed before adding split variables
   (C0) and after adding them but before the pFBA objective (C1).

7. **Stage D** — build another fresh community, minimize `sum(y_j)`, and require
   every selected product to retain at least `beta` times its Stage-C reference
   flux. Alternative globally minimum communities are enumerated with no-good
   cuts.

## Default experimental parameters

- minimum growth: `0.01`
- Stage-A production fraction `alpha`: `0.20`
- lexicographic epsilon: `1e-5`
- Stage-D product retention `beta`: `0.90`
- integrality tolerance: `1e-9`
- feasibility/optimality tolerance: `1e-6`

## Numerical status

Earlier implementations reused the Stage-A/B MILP in Stage C and sometimes
reported pFBA infeasibility even though the pre-pFBA problem was feasible.
Because adding the exact split representation `v = v+ - v-` should not alter
feasibility, this may be numerical. The current fresh-Stage-C implementation is
a diagnostic refactor and **has not yet been validated on Leftraru** at this
checkpoint.

The last validated lexicographic formulation before the C/N/S/P filter solved
successfully. After introducing the C/N/S/P filter, the inherited-solver Stage C
again reported infeasibility. The fresh Stage-C C0/C1/C2 checks are intended to
localize this behavior.

## Main files

- `misosoup/library/product_reference.py`
- `misosoup/library/product_filter.py`
- `misosoup/library/product_selection.py`
- `misosoup/library/minimal_product_communities.py`
- `test_simultaneous_products.py`
- `simultaneous_products_test.slurm`

## Running on Leftraru

```bash
cd ~/misosoup_product
python -m py_compile misosoup/library/product_reference.py
python -m py_compile misosoup/library/product_filter.py
python -m py_compile misosoup/library/product_selection.py
python -m py_compile misosoup/library/minimal_product_communities.py
python -m py_compile test_simultaneous_products.py
sbatch simultaneous_products_test.slurm
```
