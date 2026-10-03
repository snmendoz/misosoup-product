#!/bin/bash
# Submit the staged 159-sample ocean workflow.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo "Usage: $0 MANIFEST OUTPUT_ROOT [MEDIUM_FILE]" >&2
  exit 2
fi

MANIFEST="$(readlink -f "$1")"
OUTPUT_ROOT="$(readlink -m "$2")"
MEDIUM_FILE="${3:-$HOME/misosoup_product/experiments/ocean_samples/media/ocean_complete_media.csv}"
MEDIUM_FILE="$(readlink -f "$MEDIUM_FILE")"

mkdir -p "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/samples"

N_SAMPLES="$(
python - "$MANIFEST" <<'PY'
import sys, yaml
with open(sys.argv[1], "r", encoding="utf-8") as h:
    m = yaml.safe_load(h)
print(len(m["samples"]))
PY
)"

ARRAY_LAST=$((N_SAMPLES - 1))
ARRAY_RANGE="${OCEAN_ARRAY_RANGE:-0-${ARRAY_LAST}}"

# Product Scan has no artificial array throttle: SLURM decides how many
# sample jobs can run at once.  Gurobi stages remain limited to the WLS
# concurrency available to this project.
PRODUCT_ARRAY_SPEC="$ARRAY_RANGE"
GUROBI_MAX_CONCURRENT="${OCEAN_GUROBI_MAX_CONCURRENT:-2}"
GUROBI_ARRAY_SPEC="${ARRAY_RANGE}%${GUROBI_MAX_CONCURRENT}"

if (( GUROBI_MAX_CONCURRENT > 2 )); then
  echo "ERROR: current Gurobi WLS baseline supports only 2 concurrent sessions." >&2
  exit 1
fi

export OCEAN_MANIFEST="$MANIFEST"
export OCEAN_OUTPUT_ROOT="$OUTPUT_ROOT"
export OCEAN_MEDIUM_FILE="$MEDIUM_FILE"
export OCEAN_ALPHA="${OCEAN_ALPHA:-0.20}"
export OCEAN_MINIMAL_GROWTH="${OCEAN_MINIMAL_GROWTH:-0.01}"
export OCEAN_PRODUCT_RETENTION="${OCEAN_PRODUCT_RETENTION:-0.90}"
export OCEAN_LEX_TOL="${OCEAN_LEX_TOL:-1e-5}"
export OCEAN_PROD_TOL="${OCEAN_PROD_TOL:-1e-6}"
export OCEAN_INT_TOL="${OCEAN_INT_TOL:-1e-9}"
export OCEAN_UPTAKE_BOUND="${OCEAN_UPTAKE_BOUND:--1000}"
export OCEAN_MAX_MIN_COMMS="${OCEAN_MAX_MIN_COMMS:-100}"
export OCEAN_KEEP_OXYGEN_FLAG="${OCEAN_KEEP_OXYGEN_FLAG:---keep-oxygen}"

MAIL_USER="${OCEAN_MAIL_USER:-smendoza@cmm.uchile.cl}"

QC_JOB="$(sbatch --parsable \
  --job-name=ocean_qc \
  --partition="${OCEAN_QC_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_QC_CPUS:-2}" \
  --mem="${OCEAN_QC_MEM:-8G}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/qc_%j.out" \
  --error="$OUTPUT_ROOT/logs/qc_%j.err" \
  "$SCRIPT_DIR/global_qc.slurm")"
QC_JOB="${QC_JOB%%;*}"

PRODUCT_JOB="$(sbatch --parsable \
  --dependency="afterok:${QC_JOB}" \
  --job-name=ocean_product_scan \
  --partition="${OCEAN_PRODUCT_SCAN_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_PRODUCT_SCAN_CPUS:-8}" \
  --mem="${OCEAN_PRODUCT_SCAN_MEM:-32G}" \
  --array="$PRODUCT_ARRAY_SPEC" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/product_scan_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/product_scan_%A_%a.err" \
  "$SCRIPT_DIR/product_scan_parallel_array.slurm")"
PRODUCT_JOB="${PRODUCT_JOB%%;*}"

AB_JOB="$(sbatch --parsable \
  --dependency="afterok:${PRODUCT_JOB}" \
  --job-name=ocean_ab \
  --partition="${OCEAN_AB_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_AB_CPUS:-2}" \
  --mem="${OCEAN_AB_MEM:-16G}" \
  --array="$GUROBI_ARRAY_SPEC" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/ab_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/ab_%A_%a.err" \
  "$SCRIPT_DIR/reference_ab_array.slurm")"
AB_JOB="${AB_JOB%%;*}"

C_JOB="$(sbatch --parsable \
  --dependency="afterany:${AB_JOB}" \
  --job-name=ocean_c \
  --partition="${OCEAN_C_PARTITION:-largemem}" \
  --cpus-per-task="${OCEAN_C_CPUS:-2}" \
  --mem="${OCEAN_C_MEM:-32G}" \
  --array="$GUROBI_ARRAY_SPEC" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/c_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/c_%A_%a.err" \
  "$SCRIPT_DIR/stage_c_array.slurm")"
C_JOB="${C_JOB%%;*}"

D_JOB="$(sbatch --parsable \
  --dependency="afterany:${C_JOB}" \
  --job-name=ocean_d \
  --partition="${OCEAN_D_PARTITION:-largemem}" \
  --cpus-per-task="${OCEAN_D_CPUS:-2}" \
  --mem="${OCEAN_D_MEM:-24G}" \
  --array="$GUROBI_ARRAY_SPEC" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/d_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/d_%A_%a.err" \
  "$SCRIPT_DIR/stage_d_array.slurm")"
D_JOB="${D_JOB%%;*}"

COLLECT_JOB="$(sbatch --parsable \
  --dependency="afterany:${D_JOB}" \
  --job-name=ocean_collect \
  --partition=debug \
  --cpus-per-task=1 \
  --mem=2G \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/collect_%j.out" \
  --error="$OUTPUT_ROOT/logs/collect_%j.err" \
  "$SCRIPT_DIR/collect_results.slurm")"
COLLECT_JOB="${COLLECT_JOB%%;*}"

cat > "$OUTPUT_ROOT/submission_info.txt" <<EOF
manifest=$MANIFEST
output_root=$OUTPUT_ROOT
medium_file=$OCEAN_MEDIUM_FILE
number_samples=$N_SAMPLES
product_array_spec=$PRODUCT_ARRAY_SPEC
gurobi_array_spec=$GUROBI_ARRAY_SPEC
qc_job=$QC_JOB
product_scan_job=$PRODUCT_JOB
reference_ab_job=$AB_JOB
stage_c_job=$C_JOB
stage_d_job=$D_JOB
collector_job=$COLLECT_JOB
qc_resources=${OCEAN_QC_CPUS:-2}CPU/${OCEAN_QC_MEM:-8G}
product_scan_resources=${OCEAN_PRODUCT_SCAN_CPUS:-8}CPU/${OCEAN_PRODUCT_SCAN_MEM:-32G}
ab_resources=${OCEAN_AB_CPUS:-2}CPU/${OCEAN_AB_MEM:-16G}
c_resources=${OCEAN_C_CPUS:-2}CPU/${OCEAN_C_MEM:-32G}
d_resources=${OCEAN_D_CPUS:-2}CPU/${OCEAN_D_MEM:-24G}
explicit_walltime=none
EOF

echo "============================================================"
echo "Staged ocean pipeline submitted"
echo "QC             : $QC_JOB"
echo "Product Scan   : $PRODUCT_JOB  (${OCEAN_PRODUCT_SCAN_CPUS:-8} workers/sample)"
echo "Stage A+B      : $AB_JOB"
echo "Stage C        : $C_JOB"
echo "Stage D        : $D_JOB"
echo "Collector      : $COLLECT_JOB"
echo "Product array  : $PRODUCT_ARRAY_SPEC (no artificial throttle)"
echo "Gurobi array   : $GUROBI_ARRAY_SPEC"
echo "Explicit time  : none"
echo "Output         : $OUTPUT_ROOT"
echo "============================================================"
echo "squeue -j $QC_JOB,$PRODUCT_JOB,$AB_JOB,$C_JOB,$D_JOB,$COLLECT_JOB"
