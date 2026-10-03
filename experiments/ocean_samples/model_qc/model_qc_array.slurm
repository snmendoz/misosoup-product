#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J ocean_qc
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o ocean_qc_%A_%a.out
#SBATCH -e ocean_qc_%A_%a.err

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

: "${OCEAN_QC_MANIFEST:?OCEAN_QC_MANIFEST is required}"
: "${OCEAN_QC_OUTPUT_ROOT:?OCEAN_QC_OUTPUT_ROOT is required}"
: "${OCEAN_QC_MEDIUM_FILE:?OCEAN_QC_MEDIUM_FILE is required}"
: "${OCEAN_QC_CHUNK_SIZE:?OCEAN_QC_CHUNK_SIZE is required}"

OCEAN_QC_MINIMAL_GROWTH="${OCEAN_QC_MINIMAL_GROWTH:-0.01}"
OCEAN_QC_UPTAKE_BOUND="${OCEAN_QC_UPTAKE_BOUND:--1000}"
OCEAN_QC_RICH_UPTAKE_BOUND="${OCEAN_QC_RICH_UPTAKE_BOUND:--1000}"
OCEAN_QC_PROD_TOL="${OCEAN_QC_PROD_TOL:-1e-6}"
OCEAN_QC_INT_TOL="${OCEAN_QC_INT_TOL:-1e-9}"

TASK_ID="${SLURM_ARRAY_TASK_ID}"
WORK_DIR="${OCEAN_QC_OUTPUT_ROOT}/work/chunk_${TASK_ID}"

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

cat > gurobi.env <<EOF
Threads ${SLURM_CPUS_PER_TASK}
EOF

echo "============================================================"
echo "Ocean individual-model QC"
echo "Array job : ${SLURM_ARRAY_JOB_ID}"
echo "Chunk     : ${TASK_ID}"
echo "Node      : $(hostname)"
echo "Start     : $(date)"
echo "Manifest  : ${OCEAN_QC_MANIFEST}"
echo "Medium    : ${OCEAN_QC_MEDIUM_FILE}"
echo "Chunk size: ${OCEAN_QC_CHUNK_SIZE}"
echo "============================================================"

srun python -u \
  "$HOME/misosoup_product/experiments/ocean_samples/model_qc/run_model_qc_chunk.py" \
  --manifest "$OCEAN_QC_MANIFEST" \
  --chunk-index "$TASK_ID" \
  --chunk-size "$OCEAN_QC_CHUNK_SIZE" \
  --output-root "$OCEAN_QC_OUTPUT_ROOT" \
  --medium-file "$OCEAN_QC_MEDIUM_FILE" \
  --minimal-growth "$OCEAN_QC_MINIMAL_GROWTH" \
  --uptake-bound "$OCEAN_QC_UPTAKE_BOUND" \
  --rich-uptake-bound "$OCEAN_QC_RICH_UPTAKE_BOUND" \
  --production-tolerance "$OCEAN_QC_PROD_TOL" \
  --integer-tolerance "$OCEAN_QC_INT_TOL"

echo "============================================================"
echo "Finished: $(date)"
echo "============================================================"
