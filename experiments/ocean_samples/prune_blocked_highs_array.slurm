#!/bin/bash
set -euo pipefail

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

: "${OCEAN_PRUNE_LIST:?OCEAN_PRUNE_LIST is required}"
: "${OCEAN_PRUNE_OUTPUT_DIR:?OCEAN_PRUNE_OUTPUT_DIR is required}"

LINE_NUMBER=$((SLURM_ARRAY_TASK_ID + 1))
SOURCE_MODEL=$(sed -n "${LINE_NUMBER}p" "$OCEAN_PRUNE_LIST")

if [[ -z "$SOURCE_MODEL" ]]; then
  echo "No source model for array index $SLURM_ARRAY_TASK_ID" >&2
  exit 2
fi

cd "$HOME/misosoup_product"

python -u experiments/ocean_samples/prune_blocked_model_highs.py \
  --source-model "$SOURCE_MODEL" \
  --output-dir "$OCEAN_PRUNE_OUTPUT_DIR" \
  --blocked-tolerance "${OCEAN_PRUNE_TOLERANCE:-1e-9}" \
  --open-exchange-bound "${OCEAN_PRUNE_OPEN_BOUND:-1000}" \
  --progress-every "${OCEAN_PRUNE_PROGRESS_EVERY:-250}"
