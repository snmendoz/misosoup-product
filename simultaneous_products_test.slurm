#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J simult_prod
#SBATCH -p debug
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500
#SBATCH --time=00:30:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o simult_prod_%j.out
#SBATCH -e simult_prod_%j.err

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

RUN_DIR="$HOME/misosoup_runs/simultaneous_products_${SLURM_JOB_ID}"

mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

cat > gurobi.env <<EOF2
Threads ${SLURM_CPUS_PER_TASK}
EOF2

echo "=============================================="
echo "Simultaneous product test"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "=============================================="

srun python -u "$HOME/misosoup_product/test_simultaneous_products.py"

echo "=============================================="
echo "Finished: $(date)"
echo "=============================================="
