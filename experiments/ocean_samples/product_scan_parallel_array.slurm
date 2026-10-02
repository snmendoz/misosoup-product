#!/bin/bash
set -euo pipefail

ml purge
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

# Each exchange solve is one HiGHS process.  Prevent nested BLAS/OpenMP
# parallelism so SLURM_CPUS_PER_TASK is the true concurrency budget.
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

cd "$HOME/misosoup_product"

python -u experiments/ocean_samples/run_product_scan_parallel.py \
  --manifest "$OCEAN_MANIFEST" \
  --sample-index "$SLURM_ARRAY_TASK_ID" \
  --output-root "$OCEAN_OUTPUT_ROOT" \
  --medium-file "$OCEAN_MEDIUM_FILE" \
  --minimal-growth "$OCEAN_MINIMAL_GROWTH" \
  --production-tolerance "$OCEAN_PROD_TOL" \
  --uptake-bound "$OCEAN_UPTAKE_BOUND" \
  --workers "$SLURM_CPUS_PER_TASK" \
  $OCEAN_KEEP_OXYGEN_FLAG
