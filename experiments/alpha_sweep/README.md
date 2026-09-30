# Alpha sensitivity sweep

This experiment runs the validated product-preservation workflow over:

```text
alpha = 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50
```

The biological and numerical settings are otherwise kept equal to the
validated reference run:

- models: `A1R12`, `I2R16`, `I3M07`
- medium: `ac`
- minimum growth: `0.01`
- C/N/S/P product filter with oxygen retained
- lexicographic tolerance: `1e-5`
- Stage-D product retention: `0.90`
- feasibility/optimality tolerance: `1e-6`
- integrality tolerance: `1e-9`

## Design

The individual product maxima `M_i` are independent of alpha, so the script
computes the full-community product scan and C/N/S/P filtering only once.

For each alpha, it then creates fresh solver objects and runs:

1. Stage A: maximize the number of simultaneous products, giving `K*(alpha)`.
2. Stage B: maximize normalized production at fixed `K*`, giving `Q*(alpha)`
   and the selected product set.
3. Stage C: rebuild the reference state in a fresh community and solve the
   validated finite-bound split-flux pFBA.
4. Stage D: find and enumerate minimum-cardinality communities that retain at
   least 90% of every Stage-C reference product flux.

This avoids repeating the approximately 45-second individual-product scan nine
times while keeping the alpha-dependent optimization problems independent.

## Outputs

The job creates a directory:

```text
~/misosoup_runs/alpha_sweep_<JOBID>/results/
```

with:

- `alpha_sweep_results.yaml`: detailed results and product identities.
- `alpha_sweep_summary.csv`: one row per alpha, intended for plotting and
  quick inspection.
- `alpha_<value>/`: per-alpha diagnostic directory. If a stage is infeasible,
  its IIS file is written here.

The CSV contains, among other fields:

- `alpha`
- `K*`
- `Q*`
- selected products
- Stage-C pFBA objective
- minimum community size `N_min`
- number of minimum communities
- minimum-community identities
- stage timings

The YAML is checkpointed after every alpha. If one alpha fails, the script
records the error, continues with the remaining alpha values, preserves all
partial results, and exits non-zero at the end.

## Run on Leftraru

From the repository root:

```bash
cd ~/misosoup_product
sbatch experiments/alpha_sweep/alpha_sweep.slurm
```

Check status:

```bash
sacct -j <JOBID> --format=JobID,JobName,Partition,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

Inspect the compact progress log:

```bash
grep -E "alpha =|Result:|ALPHA SWEEP COMPLETE|Successful:|Failed" alpha_sweep_<JOBID>.out
```

After completion:

```bash
column -s, -t < ~/misosoup_runs/alpha_sweep_<JOBID>/results/alpha_sweep_summary.csv
```

## Custom alpha values

The Python script also accepts an explicit list:

```bash
python experiments/alpha_sweep/run_alpha_sweep.py \
  --alphas 0.05 0.10 0.20 0.30 \
  --output-dir /path/to/results
```
