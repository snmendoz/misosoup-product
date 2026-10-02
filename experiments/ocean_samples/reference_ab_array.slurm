#!/bin/bash
set -euo pipefail
ml purge
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup
export GRB_LICENSE_FILE="$HOME/.gurobi/gurobi.lic"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
cd "$HOME/misosoup_product"

python -u experiments/ocean_samples/run_reference_ab.py \
  --manifest "$OCEAN_MANIFEST" \
  --sample-index "$SLURM_ARRAY_TASK_ID" \
  --output-root "$OCEAN_OUTPUT_ROOT" \
  --medium-file "$OCEAN_MEDIUM_FILE" \
  --alpha "$OCEAN_ALPHA" \
  --minimal-growth "$OCEAN_MINIMAL_GROWTH" \
  --production-tolerance "$OCEAN_PROD_TOL" \
  --integer-tolerance "$OCEAN_INT_TOL" \
  --uptake-bound "$OCEAN_UPTAKE_BOUND" \
  $OCEAN_KEEP_OXYGEN_FLAG
