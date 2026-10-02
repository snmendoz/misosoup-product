#!/bin/bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: $0 MANIFEST OUTPUT_ROOT MEDIUM_FILE" >&2
  exit 2
fi

MANIFEST="$(readlink -f "$1")"
OUTPUT_ROOT="$(readlink -m "$2")"
MEDIUM_FILE="$(readlink -f "$3")"

mkdir -p "$OUTPUT_ROOT/logs" "$OUTPUT_ROOT/samples"

export OCEAN_MANIFEST="$MANIFEST"
export OCEAN_OUTPUT_ROOT="$OUTPUT_ROOT"
export OCEAN_MEDIUM_FILE="$MEDIUM_FILE"
export OCEAN_MINIMAL_GROWTH="${OCEAN_MINIMAL_GROWTH:-0.01}"
export OCEAN_PROD_TOL="${OCEAN_PROD_TOL:-1e-6}"
export OCEAN_UPTAKE_BOUND="${OCEAN_UPTAKE_BOUND:--1000}"
export OCEAN_KEEP_OXYGEN_FLAG="${OCEAN_KEEP_OXYGEN_FLAG:---keep-oxygen}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAIL_USER="${OCEAN_MAIL_USER:-smendoza@cmm.uchile.cl}"

LAST_INDEX="$(python - "$MANIFEST" <<'PY'
import sys
import yaml
with open(sys.argv[1], encoding="utf-8") as handle:
    payload = yaml.safe_load(handle)
print(len(payload["samples"]) - 1)
PY
)"

JOB_ID="$(sbatch --parsable \
  --job-name=ocean_product_scan \
  --partition="${OCEAN_PRODUCT_SCAN_PARTITION:-general}" \
  --cpus-per-task="${OCEAN_PRODUCT_SCAN_CPUS:-8}" \
  --mem="${OCEAN_PRODUCT_SCAN_MEM:-32G}" \
  --array="0-${LAST_INDEX}" \
  --mail-user="$MAIL_USER" --mail-type=FAIL,END \
  --output="$OUTPUT_ROOT/logs/product_scan_%A_%a.out" \
  --error="$OUTPUT_ROOT/logs/product_scan_%A_%a.err" \
  "$SCRIPT_DIR/product_scan_parallel_array.slurm")"
JOB_ID="${JOB_ID%%;*}"

cat > "$OUTPUT_ROOT/product_scan_submission_info.txt" <<EOF
manifest=$MANIFEST
output_root=$OUTPUT_ROOT
medium_file=$MEDIUM_FILE
product_scan_job=$JOB_ID
sample_array=0-$LAST_INDEX
cpus_per_sample=${OCEAN_PRODUCT_SCAN_CPUS:-8}
memory_per_sample=${OCEAN_PRODUCT_SCAN_MEM:-32G}
solver=HiGHS
exchange_parallelism=SLURM_CPUS_PER_TASK
explicit_walltime=none
EOF

echo "Product Scan array : $JOB_ID"
echo "Samples            : 0-$LAST_INDEX"
echo "CPUs/sample        : ${OCEAN_PRODUCT_SCAN_CPUS:-8}"
echo "Memory/sample      : ${OCEAN_PRODUCT_SCAN_MEM:-32G}"
echo "Exchange workers   : ${OCEAN_PRODUCT_SCAN_CPUS:-8} per running sample"
echo "Explicit time limit: none"
echo "Output             : $OUTPUT_ROOT"
