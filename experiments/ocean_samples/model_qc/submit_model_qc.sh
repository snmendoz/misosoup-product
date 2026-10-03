#!/bin/bash
# Run individual growth QC for all 1375 ocean MAG models or a named subset.

set -e
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OCEAN_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

DATA_REPO="${OCEAN_DATA_REPO:-$HOME/tara_chile_metabolic_models}"
DATA_REPO="${DATA_REPO/#\~/$HOME}"

MATRIX="${OCEAN_ABUNDANCE_MATRIX:-$DATA_REPO/Data/CEODOS_MAG_TPMs_IDs_matrix.tsv}"
MODELS_DIR="${OCEAN_MODELS_DIR:-$DATA_REPO/Models/gapseq/models}"

TARA_MEDIUM="$DATA_REPO/Models/gapseq/media/mathomics.txt"
BUNDLED_MEDIUM="$OCEAN_DIR/media/complex_media_gapseq2.csv"

if [[ -n "${OCEAN_MEDIUM_FILE:-}" ]]; then
    MEDIUM_FILE="$OCEAN_MEDIUM_FILE"
elif [[ -f "$TARA_MEDIUM" ]]; then
    MEDIUM_FILE="$TARA_MEDIUM"
else
    MEDIUM_FILE="$BUNDLED_MEDIUM"
fi

MAG_IDS="${OCEAN_QC_MAG_IDS:-}"
CHUNK_SIZE="${OCEAN_QC_CHUNK_SIZE:-25}"
MINIMAL_GROWTH="${OCEAN_QC_MINIMAL_GROWTH:-0.01}"
UPTAKE_BOUND="${OCEAN_QC_UPTAKE_BOUND:--1000}"
RICH_UPTAKE_BOUND="${OCEAN_QC_RICH_UPTAKE_BOUND:--1000}"

PARTITION="${OCEAN_QC_PARTITION:-general}"
TIME_LIMIT="${OCEAN_QC_TIME_LIMIT:-01:00:00}"
CPUS_PER_TASK="${OCEAN_QC_CPUS_PER_TASK:-2}"
MEMORY="${OCEAN_QC_MEMORY:-4G}"
MAX_CONCURRENT="${OCEAN_QC_MAX_CONCURRENT:-2}"
WLS_SESSIONS="${OCEAN_WLS_SESSIONS:-2}"

if (( CHUNK_SIZE < 1 )); then
    echo "ERROR: OCEAN_QC_CHUNK_SIZE must be >= 1" >&2
    exit 1
fi

if (( MAX_CONCURRENT > WLS_SESSIONS )); then
    echo "ERROR: OCEAN_QC_MAX_CONCURRENT=$MAX_CONCURRENT exceeds" >&2
    echo "configured WLS sessions=$WLS_SESSIONS." >&2
    exit 1
fi

if [[ ! -f "$MATRIX" ]]; then
    echo "ERROR: abundance matrix not found: $MATRIX" >&2
    exit 1
fi

if [[ ! -d "$MODELS_DIR" ]]; then
    echo "ERROR: models directory not found: $MODELS_DIR" >&2
    exit 1
fi

if [[ ! -f "$MEDIUM_FILE" ]]; then
    echo "ERROR: medium file not found: $MEDIUM_FILE" >&2
    exit 1
fi

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="${OCEAN_QC_RUN_ROOT:-$HOME/misosoup_runs/ocean_model_qc_${RUN_STAMP}}"
MANIFEST="$RUN_ROOT/model_qc_manifest.yaml"

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/models" "$RUN_ROOT/chunks" "$RUN_ROOT/work"

PREPARE_ARGS=(
    --matrix "$MATRIX"
    --models-dir "$MODELS_DIR"
    --output "$MANIFEST"
    --expected-samples 159
    --expected-mags 1375
)

if [[ -n "$MAG_IDS" ]]; then
    PREPARE_ARGS+=(--mag-ids "$MAG_IDS")
fi

echo "============================================================"
echo "Preparing individual-model QC"
echo "Matrix             : $MATRIX"
echo "Models             : $MODELS_DIR"
echo "Medium             : $MEDIUM_FILE"
echo "MAG subset         : ${MAG_IDS:-ALL}"
echo "Minimal growth     : $MINIMAL_GROWTH"
echo "Rich uptake bound  : $RICH_UPTAKE_BOUND"
echo "Chunk size         : $CHUNK_SIZE"
echo "Max concurrent     : $MAX_CONCURRENT"
echo "Run root           : $RUN_ROOT"
echo "============================================================"

python -u "$SCRIPT_DIR/prepare_model_qc_manifest.py" "${PREPARE_ARGS[@]}"

N_MODELS="$(
python - "$MANIFEST" <<'PY'
import sys
import yaml
with open(sys.argv[1], encoding="utf-8") as handle:
    manifest = yaml.safe_load(handle)
print(manifest["number_models"])
PY
)"

N_CHUNKS="$(( (N_MODELS + CHUNK_SIZE - 1) / CHUNK_SIZE ))"
LAST_CHUNK="$(( N_CHUNKS - 1 ))"
ARRAY_RANGE="${OCEAN_QC_ARRAY_RANGE:-0-${LAST_CHUNK}}"
ARRAY_SPEC="${ARRAY_RANGE}%${MAX_CONCURRENT}"

export OCEAN_QC_MANIFEST="$MANIFEST"
export OCEAN_QC_OUTPUT_ROOT="$RUN_ROOT"
export OCEAN_QC_MEDIUM_FILE="$MEDIUM_FILE"
export OCEAN_QC_CHUNK_SIZE="$CHUNK_SIZE"
export OCEAN_QC_MINIMAL_GROWTH="$MINIMAL_GROWTH"
export OCEAN_QC_UPTAKE_BOUND="$UPTAKE_BOUND"
export OCEAN_QC_RICH_UPTAKE_BOUND="$RICH_UPTAKE_BOUND"

ARRAY_JOB_ID="$(
sbatch --parsable \
  --export=ALL \
  --partition="$PARTITION" \
  --time="$TIME_LIMIT" \
  --cpus-per-task="$CPUS_PER_TASK" \
  --mem="$MEMORY" \
  --array="$ARRAY_SPEC" \
  --output="$RUN_ROOT/logs/model_qc_%A_%a.out" \
  --error="$RUN_ROOT/logs/model_qc_%A_%a.err" \
  "$SCRIPT_DIR/model_qc_array.slurm"
)"

COLLECTOR_JOB_ID="$(
sbatch --parsable \
  --export=ALL \
  --dependency="afterany:${ARRAY_JOB_ID}" \
  --output="$RUN_ROOT/logs/model_qc_collect_%j.out" \
  --error="$RUN_ROOT/logs/model_qc_collect_%j.err" \
  "$SCRIPT_DIR/model_qc_collect.slurm"
)"

cat > "$RUN_ROOT/submission_info.txt" <<EOF
run_root=$RUN_ROOT
manifest=$MANIFEST
medium=$MEDIUM_FILE
number_models=$N_MODELS
chunk_size=$CHUNK_SIZE
number_chunks=$N_CHUNKS
array_spec=$ARRAY_SPEC
array_job_id=$ARRAY_JOB_ID
collector_job_id=$COLLECTOR_JOB_ID
minimal_growth=$MINIMAL_GROWTH
rich_uptake_bound=$RICH_UPTAKE_BOUND
mag_subset=$MAG_IDS
EOF

echo
echo "============================================================"
echo "Submitted individual-model QC"
echo "============================================================"
echo "Models        : $N_MODELS"
echo "Chunks        : $N_CHUNKS"
echo "Array         : $ARRAY_SPEC"
echo "Array job     : $ARRAY_JOB_ID"
echo "Collector job : $COLLECTOR_JOB_ID"
echo "Run root      : $RUN_ROOT"
echo
echo "Status:"
echo "  squeue -j $ARRAY_JOB_ID,$COLLECTOR_JOB_ID"
echo
echo "After collector:"
echo "  column -s, -t < \"$RUN_ROOT/model_qc.csv\" | less -S"
echo "  cat \"$RUN_ROOT/aggregate.yaml\""
echo "============================================================"
