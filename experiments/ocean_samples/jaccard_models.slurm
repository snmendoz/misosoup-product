#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J tara_jaccard
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o tara_jaccard_%j.out
#SBATCH -e tara_jaccard_%j.err

set -e
set -u
set -o pipefail

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

JACCARD_DATA_ROOT="${JACCARD_DATA_ROOT:-$HOME/tara_chile_metabolic_models}"
: "${JACCARD_OUTPUT_ROOT:?JACCARD_OUTPUT_ROOT is required}"

echo "============================================================"
echo "TARA Chile Jaccard matrices"
echo "Job       : ${SLURM_JOB_ID}"
echo "Node      : $(hostname)"
echo "Start     : $(date)"
echo "Data root : ${JACCARD_DATA_ROOT}"
echo "Output    : ${JACCARD_OUTPUT_ROOT}"
echo "============================================================"

srun python -u   "$HOME/misosoup_product/experiments/ocean_samples/jaccard_models.py"   --data-root "$JACCARD_DATA_ROOT"   --output-root "$JACCARD_OUTPUT_ROOT"   --expected-models 1375   --progress-every 100

echo "============================================================"
echo "Finished: $(date)"
echo "Results : ${JACCARD_OUTPUT_ROOT}"
echo "============================================================"
