#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J ocean_sample
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o ocean_sample_%A_%a.out
#SBATCH -e ocean_sample_%A_%a.err

set -e
set -u
set -o pipefail

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

export GRB_LICENSE_FILE="$HOME/.gurobi/gurobi.lic"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

: "${OCEAN_MANIFEST:?OCEAN_MANIFEST is required}"
: "${OCEAN_OUTPUT_ROOT:?OCEAN_OUTPUT_ROOT is required}"
: "${OCEAN_MEDIUM_FILE:?OCEAN_MEDIUM_FILE is required}"

OCEAN_ALPHA="${OCEAN_ALPHA:-0.20}"
OCEAN_MINIMAL_GROWTH="${OCEAN_MINIMAL_GROWTH:-0.01}"
OCEAN_PRODUCT_RETENTION="${OCEAN_PRODUCT_RETENTION:-0.90}"
OCEAN_LEX_TOL="${OCEAN_LEX_TOL:-1e-5}"
OCEAN_PROD_TOL="${OCEAN_PROD_TOL:-1e-6}"
OCEAN_INT_TOL="${OCEAN_INT_TOL:-1e-9}"
OCEAN_UPTAKE_BOUND="${OCEAN_UPTAKE_BOUND:--1000}"
OCEAN_MAX_MIN_COMMS="${OCEAN_MAX_MIN_COMMS:-100}"

TASK_ID="${SLURM_ARRAY_TASK_ID}"
WORK_DIR="${OCEAN_OUTPUT_ROOT}/work/task_${TASK_ID}"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

cat > gurobi.env <<EOF2
Threads ${SLURM_CPUS_PER_TASK}
EOF2

echo "============================================================"
echo "Ocean sample array task"
echo "Array job : ${SLURM_ARRAY_JOB_ID}"
echo "Task      : ${TASK_ID}"
echo "Node      : $(hostname)"
echo "Start     : $(date)"
echo "Manifest  : ${OCEAN_MANIFEST}"
echo "Medium    : ${OCEAN_MEDIUM_FILE}"
echo "Output    : ${OCEAN_OUTPUT_ROOT}"
echo "============================================================"

srun python -u \
  "$HOME/misosoup_product/experiments/ocean_samples/run_sample.py" \
  --manifest "$OCEAN_MANIFEST" \
  --sample-index "$TASK_ID" \
  --output-root "$OCEAN_OUTPUT_ROOT" \
  --medium-file "$OCEAN_MEDIUM_FILE" \
  --alpha "$OCEAN_ALPHA" \
  --minimal-growth "$OCEAN_MINIMAL_GROWTH" \
  --product-retention "$OCEAN_PRODUCT_RETENTION" \
  --lexicographic-tolerance "$OCEAN_LEX_TOL" \
  --production-tolerance "$OCEAN_PROD_TOL" \
  --integer-tolerance "$OCEAN_INT_TOL" \
  --uptake-bound "$OCEAN_UPTAKE_BOUND" \
  --max-minimal-communities "$OCEAN_MAX_MIN_COMMS"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
