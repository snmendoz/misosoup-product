#!/bin/bash
# One-time structural model reduction: FVA with all environmental exchanges open.
set -euo pipefail

ml purge
source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

export GRB_LICENSE_FILE="$HOME/.gurobi/gurobi.lic"
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

: "${OCEAN_PRUNE_MANIFEST:?OCEAN_PRUNE_MANIFEST is required}"
: "${OCEAN_PRUNE_MODEL_DIR:?OCEAN_PRUNE_MODEL_DIR is required}"
: "${OCEAN_PRUNE_OUTPUT_MANIFEST:?OCEAN_PRUNE_OUTPUT_MANIFEST is required}"

cd "$HOME/misosoup_product"

python -u experiments/ocean_samples/prune_blocked_models.py \
  --manifest "$OCEAN_PRUNE_MANIFEST" \
  --output-model-dir "$OCEAN_PRUNE_MODEL_DIR" \
  --output-manifest "$OCEAN_PRUNE_OUTPUT_MANIFEST" \
  --summary-csv "${OCEAN_PRUNE_SUMMARY_CSV:-$OCEAN_PRUNE_MODEL_DIR/blocked_reaction_pruning_summary.csv}" \
  --workers "${OCEAN_PRUNE_WORKERS:-2}" \
  --blocked-tolerance "${OCEAN_PRUNE_TOLERANCE:-1e-9}" \
  --open-exchange-bound "${OCEAN_PRUNE_OPEN_BOUND:-1000}" \
  ${OCEAN_PRUNE_EXTRA_ARGS:-}
