#!/bin/bash
set -euo pipefail

if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: $0 SOURCE_DIR OUTPUT_DIR RUN_DIR [MAX_CONCURRENT]" >&2
  exit 2
fi

SOURCE_DIR=$(readlink -f "$1")
OUTPUT_DIR=$(readlink -m "$2")
RUN_DIR=$(readlink -m "$3")
MAX_CONCURRENT="${4:-48}"

mkdir -p "$OUTPUT_DIR" "$RUN_DIR/logs"

cd "$HOME/misosoup_product"

PYTHON="${OCEAN_PRUNE_PYTHON:-$HOME/anaconda3/envs/misosoup/bin/python}"
if [[ ! -x "$PYTHON" ]]; then
  echo "ERROR: Python executable not found: $PYTHON" >&2
  exit 1
fi

"$PYTHON" -c 'import reframed, scipy' || {
  echo "ERROR: misosoup Python environment is missing reframed or scipy." >&2
  exit 1
}

LIST_FILE="$RUN_DIR/pending_models.txt"

"$PYTHON" experiments/ocean_samples/prepare_blocked_prune_pending.py \
  --source-dir "$SOURCE_DIR" \
  --output-dir "$OUTPUT_DIR" \
  --output-list "$LIST_FILE" \
  --blocked-tolerance "${OCEAN_PRUNE_TOLERANCE:-1e-9}" \
  --open-exchange-bound "${OCEAN_PRUNE_OPEN_BOUND:-1000}"

N_PENDING=$(grep -cve '^$' "$LIST_FILE" || true)

if [[ "$N_PENDING" -eq 0 ]]; then
  echo "Nothing to submit: every current source model already has a valid reduced model."
  exit 0
fi

ARRAY_MAX=$((N_PENDING - 1))
JOB_NAME="${OCEAN_PRUNE_JOB_NAME:-prune_blocked}"

JOB_ID=$(sbatch --parsable \
  -J "$JOB_NAME" \
  -p "${OCEAN_PRUNE_PARTITION:-general}" \
  -c 1 \
  --mem="${OCEAN_PRUNE_MEMORY:-3G}" \
  --time="${OCEAN_PRUNE_TIME:-08:00:00}" \
  --array="0-${ARRAY_MAX}%${MAX_CONCURRENT}" \
  --mail-user="${OCEAN_PRUNE_MAIL_USER:-smendoza@cmm.uchile.cl}" \
  --mail-type=FAIL,END \
  -o "$RUN_DIR/logs/%A_%a.out" \
  -e "$RUN_DIR/logs/%A_%a.err" \
  --export=ALL,OCEAN_PRUNE_LIST="$LIST_FILE",OCEAN_PRUNE_OUTPUT_DIR="$OUTPUT_DIR" \
  experiments/ocean_samples/prune_blocked_highs_array.slurm)

printf "%s\n" "$JOB_ID" > "$RUN_DIR/job_id.txt"
printf "%s\n" "$SOURCE_DIR" > "$RUN_DIR/source_dir.txt"
printf "%s\n" "$OUTPUT_DIR" > "$RUN_DIR/output_dir.txt"

echo "Submitted $N_PENDING models as job array $JOB_ID"
echo "Max concurrent tasks: $MAX_CONCURRENT"
echo "Pending list: $LIST_FILE"
echo "Output models: $OUTPUT_DIR"
echo "Logs: $RUN_DIR/logs"
