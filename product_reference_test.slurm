#!/bin/bash
#--------------- Script SBATCH - NLHPC ----------------
#SBATCH -J product_ref
#SBATCH -p debug
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=1500
#SBATCH --time=00:30:00
#SBATCH --mail-user=smendoza@cmm.uchile.cl
#SBATCH --mail-type=ALL
#SBATCH -o product_ref_%j.out
#SBATCH -e product_ref_%j.err

#---------------- Toolchain ----------------------------
set -e

ml purge

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate misosoup

set -u
set -o pipefail

#---------------- Gurobi -------------------------------
export GRB_LICENSE_FILE="$HOME/.gurobi/gurobi.lic"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

#---------------- Execution ----------------------------

RUN_DIR="$HOME/misosoup_runs/product_reference_${SLURM_JOB_ID}"
mkdir -p "$RUN_DIR"
cd "$RUN_DIR"

cat > gurobi.env <<EOF2
Threads ${SLURM_CPUS_PER_TASK}
EOF2

echo "=============================================="
echo "MiSoSoup product-reference test"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Start: $(date)"
echo "=============================================="

srun python "$HOME/misosoup_product/test_product_reference.py"

echo "=============================================="
echo "Finished: $(date)"
echo "=============================================="
