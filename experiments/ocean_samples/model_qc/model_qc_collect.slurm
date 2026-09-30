#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J ocean_qc_collect
#SBATCH -p debug
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL

set -e
set -u
set -o pipefail

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

: "${OCEAN_QC_MANIFEST:?OCEAN_QC_MANIFEST is required}"
: "${OCEAN_QC_OUTPUT_ROOT:?OCEAN_QC_OUTPUT_ROOT is required}"

srun python -u \
  "$HOME/misosoup_product/experiments/ocean_samples/model_qc/collect_model_qc.py" \
  --manifest "$OCEAN_QC_MANIFEST" \
  --output-root "$OCEAN_QC_OUTPUT_ROOT"
