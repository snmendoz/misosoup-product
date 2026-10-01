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
   (C0) and after adding them but before the pFBA objective (C1). The final
   pFBA solve is C2.

   The split variables use **finite, reaction-aware bounds** derived from the
   original reaction bounds:

   `0 <= v_r+ <= max(0, ub_r)`

   `0 <= v_r- <= max(0, -lb_r)`

   Infinite original bounds are replaced using MiSoSoup's practical
   `BOUND_INF = 1000` convention.

7. **Stage D** — build another fresh community, minimize `sum(y_j)`, and require
   every selected product to retain at least `beta` times its Stage-C reference
   flux. Alternative globally minimum communities are enumerated with no-good
   cuts.

## Stage-D audit outputs

The ocean-sample workflow writes two explicit audit tables for every successful
Stage-D run.

1. `product_preservation_audit.csv` contains one row for every
   minimum-community/product pair. It records the full-community product
   maximum, the Stage-C reference flux, the Stage-D required lower bound
   (`beta * v_i^ref`), the observed Stage-D flux, the retention fraction,
   the numerical margin to the requirement, and a PASS/FAIL flag. The workflow
   raises an error if any returned minimum community violates a product
   requirement beyond the configured production tolerance.

2. `medium_uptake_audit.csv` contains one row for every
   minimum-community/medium-exchange pair. It records the configured uptake
   lower bound, the raw Stage-D exchange flux, the actual uptake magnitude
   (negative exchange flux converted to a positive magnitude), secretion when
   present, the fraction of the allowed uptake that is used, the distance from
   the lower bound, and whether the returned Stage-D solution is at the uptake
   limit.

The medium audit reports the fluxes of the particular optimal Stage-D solution
returned by the solver. Because Stage D minimizes community cardinality rather
than nutrient uptake, these fluxes should be interpreted as an audit of the
returned solution, not as proof that a nutrient uptake is uniquely required.

## Default experimental parameters

- minimum growth: `0.01`
- Stage-A production fraction `alpha`: `0.20`
- lexicographic epsilon: `1e-5`
- Stage-D product retention `beta`: `0.90`
- integrality tolerance: `1e-9`
- feasibility/optimality tolerance: `1e-6`

## Numerical status

The fresh Stage-C implementation has now been validated on Leftraru.

Earlier versions used unbounded split variables,

`0 <= v_r+ < infinity`

`0 <= v_r- < infinity`,

and showed a characteristic failure in which C0 and C1 were feasible, but C2
was reported as infeasible after introducing the pFBA objective. Because C1 and
C2 have the same feasible region and differ only in the objective, this behavior
was inconsistent with exact mathematical feasibility and pointed to a numerical
issue in the split-flux formulation.

Replacing the infinite split-variable upper bounds with finite reaction-aware
bounds resolved the C2 failure without changing Stage A, Stage B, the biological
constraints, `alpha`, `epsilon`, or solver tolerances.

The validated Leftraru run was job `13702509`.

## Validated reference run

Configuration:

- organisms: `A1R12`, `I2R16`, `I3M07`
- medium: `ac`
- minimum growth: `0.01`
- producible global exchanges before filtering: `79`
- candidate products after the C/N/S/P filter: `73`
- `alpha = 0.20`
- `epsilon = 1e-5`
- `beta = 0.90`

Stage A:

- `K* = 7`
- selected products:
  - `R_EX_5mtr_e`
  - `R_EX_co2_e`
  - `R_EX_gua_e`
  - `R_EX_h2s_e`
  - `R_EX_oxa_e`
  - `R_EX_ptrc_e`
  - `R_EX_s_e`

Stage B:

- `Q* = 1.8047669375271047`
- reconstructed `sum(q_i) = 1.8047669375271045`
- reconstruction difference: `2.22e-16`
- selected products:
  - `R_EX_5mtr_e`
  - `R_EX_co2_e`
  - `R_EX_gua_e`
  - `R_EX_h2s_e`
  - `R_EX_oxa_e`
  - `R_EX_s_e`
  - `R_EX_spmd_e`

Thus Stage B preserves the Stage-A cardinality but replaces putrescine
(`R_EX_ptrc_e`) with spermidine (`R_EX_spmd_e`) when maximizing normalized
production.

Stage C:

- C0: feasible
- C1: feasible
- C2 pFBA: optimal
- lexicographic floor: `1.8047569375271046`
- achieved `sum(q_i) = 1.8047569375271049`
- pFBA objective: `1451.2348882200026`

Stage-C reference product fluxes:

| Product | Reference flux |
|---|---:|
| `R_EX_5mtr_e` | 0.3884547604309968 |
| `R_EX_co2_e` | 6.207117806575437 |
| `R_EX_gua_e` | 0.750462422309124 |
| `R_EX_h2s_e` | 2.7715519184719253 |
| `R_EX_oxa_e` | 1.8761532124918106 |
| `R_EX_s_e` | 2.4957626988478525 |
| `R_EX_spmd_e` | 0.3884547604309968 |

Stage D:

- global minimum community size: `2`
- product retention: `0.90`
- number of minimum communities found: `1`
- enumeration reported complete
- minimum community: `I2R16 + I3M07`
- growth:
  - `I2R16 = 0.009999999999763531`
  - `I3M07 = 0.009999999999763531`
  - community growth = `0.019999999999527063`

The minimum community satisfies all seven product-retention constraints.

## Timing of validated run

- individual product scan: `45.31 s`
- Stage A: `0.37 s`
- Stage B.1: `0.37 s`
- Stage B.2: `<0.001 s`
- Stage C: `0.77 s`
- Stage D: `0.42 s`
- total: `49.19 s`

The individual-product scan dominates the runtime for this three-organism test
case.

## Next validation step

The next planned experiment is a sensitivity analysis over:

`alpha = 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50`.

For each `alpha`, record at least:

- `K*(alpha)`
- `Q*(alpha)`
- Stage-B selected product set
- Stage-C pFBA objective
- minimum product-preserving community size
- enumerated minimum communities

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
