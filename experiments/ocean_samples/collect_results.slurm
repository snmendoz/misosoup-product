#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J ocean_collect
#SBATCH -p debug
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o ocean_collect_%j.out
#SBATCH -e ocean_collect_%j.err

set -e
set -u
set -o pipefail

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

: "${OCEAN_MANIFEST:?OCEAN_MANIFEST is required}"
: "${OCEAN_OUTPUT_ROOT:?OCEAN_OUTPUT_ROOT is required}"

echo "============================================================"
echo "Collecting ocean-sample results"
echo "Job ID   : ${SLURM_JOB_ID}"
echo "Node     : $(hostname)"
echo "Start    : $(date)"
echo "Manifest : ${OCEAN_MANIFEST}"
echo "Output   : ${OCEAN_OUTPUT_ROOT}"
echo "============================================================"

srun python -u \
  "$HOME/misosoup_product/experiments/ocean_samples/collect_results.py" \
  --manifest "$OCEAN_MANIFEST" \
  --output-root "$OCEAN_OUTPUT_ROOT"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
