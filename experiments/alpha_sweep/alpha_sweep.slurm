#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J alpha_sweep
#SBATCH -p debug
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500
#SBATCH --time=00:30:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o alpha_sweep_%j.out
#SBATCH -e alpha_sweep_%j.err

set -e

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

set -u
set -o pipefail

export GRB_LICENSE_FILE="$HOME/.gurobi/gurobi.lic"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

RUN_DIR="$HOME/misosoup_runs/alpha_sweep_${SLURM_JOB_ID}"
RESULT_DIR="${RUN_DIR}/results"

mkdir -p "${RESULT_DIR}"
cd "${RUN_DIR}"

cat > gurobi.env <<EOF2
Threads ${SLURM_CPUS_PER_TASK}
EOF2

echo "=============================================="
echo "MiSoSoup alpha sensitivity sweep"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "Results: ${RESULT_DIR}"
echo "=============================================="

srun python -u \
  "$HOME/misosoup_product/experiments/alpha_sweep/run_alpha_sweep.py" \
  --output-dir "${RESULT_DIR}"

echo "=============================================="
echo "Finished: $(date)"
echo "YAML: ${RESULT_DIR}/alpha_sweep_results.yaml"
echo "CSV : ${RESULT_DIR}/alpha_sweep_summary.csv"
echo "=============================================="
