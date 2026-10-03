#!/bin/bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 MANIFEST SAMPLE_INDEX OUTPUT_ROOT MEDIUM_FILE" >&2
  exit 2
fi

MANIFEST="$(readlink -f "$1")"
SAMPLE_INDEX="$2"
OUTPUT_ROOT="$(readlink -m "$3")"
MEDIUM_FILE="$(readlink -f "$4")"

mkdir -p "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/samples"

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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIL_USER="${OCEAN_MAIL_USER:-smendoza@cmm.uchile.cl}"

PRODUCT_JOB="$(sbatch --parsable \
  --job-name=ocean_ps_s${SAMPLE_INDEX} \
  --partition="${OCEAN_PRODUCT_SCAN_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_PRODUCT_SCAN_CPUS:-8}" \
  --mem="${OCEAN_PRODUCT_SCAN_MEM:-32G}" \
  --array="${SAMPLE_INDEX}-${SAMPLE_INDEX}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/product_scan_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/product_scan_%A_%a.err" \
  "$SCRIPT_DIR/product_scan_parallel_array.slurm")"
PRODUCT_JOB="${PRODUCT_JOB%%;*}"

AB_JOB="$(sbatch --parsable \
  --dependency="afterok:${PRODUCT_JOB}" \
  --job-name=ocean_ab_s${SAMPLE_INDEX} \
  --partition="${OCEAN_AB_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_AB_CPUS:-2}" \
  --mem="${OCEAN_AB_MEM:-16G}" \
  --array="${SAMPLE_INDEX}-${SAMPLE_INDEX}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/ab_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/ab_%A_%a.err" \
  "$SCRIPT_DIR/reference_ab_array.slurm")"
AB_JOB="${AB_JOB%%;*}"

C_JOB="$(sbatch --parsable \
  --dependency="afterok:${AB_JOB}" \
  --job-name=ocean_c_s${SAMPLE_INDEX} \
  --partition="${OCEAN_C_PARTITION:-largemem}" \
  --cpus-per-task="${OCEAN_C_CPUS:-2}" \
  --mem="${OCEAN_C_MEM:-32G}" \
  --array="${SAMPLE_INDEX}-${SAMPLE_INDEX}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/c_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/c_%A_%a.err" \
  "$SCRIPT_DIR/stage_c_array.slurm")"
C_JOB="${C_JOB%%;*}"

D_JOB="$(sbatch --parsable \
  --dependency="afterok:${C_JOB}" \
  --job-name=ocean_d_s${SAMPLE_INDEX} \
  --partition="${OCEAN_D_PARTITION:-largemem}" \
  --cpus-per-task="${OCEAN_D_CPUS:-2}" \
  --mem="${OCEAN_D_MEM:-24G}" \
  --array="${SAMPLE_INDEX}-${SAMPLE_INDEX}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/d_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/d_%A_%a.err" \
  "$SCRIPT_DIR/stage_d_array.slurm")"
D_JOB="${D_JOB%%;*}"

cat > "$OUTPUT_ROOT/submission_info.txt" <<EOF
manifest=$MANIFEST
sample_index=$SAMPLE_INDEX
output_root=$OUTPUT_ROOT
medium_file=$MEDIUM_FILE
product_scan_job=$PRODUCT_JOB
ab_job=$AB_JOB
stage_c_job=$C_JOB
stage_d_job=$D_JOB
product_scan_resources=${OCEAN_PRODUCT_SCAN_CPUS:-8}CPU/${OCEAN_PRODUCT_SCAN_MEM:-32G}
ab_resources=${OCEAN_AB_CPUS:-2}CPU/${OCEAN_AB_MEM:-16G}
c_resources=${OCEAN_C_CPUS:-2}CPU/${OCEAN_C_MEM:-32G}
d_resources=${OCEAN_D_CPUS:-2}CPU/${OCEAN_D_MEM:-24G}
explicit_walltime=none
EOF

echo "Product Scan      : $PRODUCT_JOB  (${OCEAN_PRODUCT_SCAN_CPUS:-8} CPU, ${OCEAN_PRODUCT_SCAN_MEM:-32G})"
echo "Stage A+B         : $AB_JOB  (${OCEAN_AB_CPUS:-2} CPU, ${OCEAN_AB_MEM:-16G})"
echo "Stage C pFBA      : $C_JOB  (${OCEAN_C_CPUS:-2} CPU, ${OCEAN_C_MEM:-32G})"
echo "Stage D min comm  : $D_JOB  (${OCEAN_D_CPUS:-2} CPU, ${OCEAN_D_MEM:-24G})"
echo "Explicit time limit: none"
echo "Output            : $OUTPUT_ROOT"
