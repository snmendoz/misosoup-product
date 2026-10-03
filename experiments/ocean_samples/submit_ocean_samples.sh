#!/bin/bash
# Submit the 159 ocean samples as a Slurm job array plus a dependent collector.

set -e
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

usage() {
    cat <<'EOF'
Usage:
  submit_ocean_samples.sh [ABUNDANCE_MATRIX] [TARA_REPO]

By default, the abundance matrix is read from:
  <TARA_REPO>/Data/CEODOS_MAG_TPMs_IDs_matrix.tsv

Example using the repository defaults:
  bash experiments/ocean_samples/submit_ocean_samples.sh

Example with explicit paths:
  bash experiments/ocean_samples/submit_ocean_samples.sh \
    /path/to/CEODOS_MAG_TPMs_IDs_matrix.tsv \
    "$HOME/tara_chile_metabolic_models"

Important environment overrides:
  OCEAN_PRESENCE_THRESHOLD=0
  OCEAN_ALPHA=0.20
  OCEAN_MINIMAL_GROWTH=0.01
  OCEAN_PRODUCT_RETENTION=0.90
  OCEAN_UPTAKE_BOUND=-1000
  OCEAN_MAX_MIN_COMMS=100
  OCEAN_MAX_MAGS_PER_SAMPLE=4
  OCEAN_ALLOW_MISSING_MODELS=1

Slurm:
  OCEAN_PARTITION=general
  OCEAN_TIME_LIMIT=12:00:00
  OCEAN_CPUS_PER_TASK=4
  OCEAN_MEMORY=32G
  OCEAN_MAX_CONCURRENT=2
  OCEAN_ARRAY_RANGE=0-158   # optional pilot, e.g. 0-3

The default concurrency is 2 because the current Gurobi WLS academic license
supports two concurrent sessions. Do not increase it unless the license limit
has changed.

Optional explicit MAG/model mapping:
  OCEAN_MODEL_MAP=/path/to/mag_model_map.tsv
EOF
}

if [[ $# -gt 2 ]]; then
    usage
    exit 2
fi

DATA_REPO="${2:-${OCEAN_DATA_REPO:-$HOME/tara_chile_metabolic_models}}"
DATA_REPO="${DATA_REPO/#\~/$HOME}"

MODELS_DIR="${OCEAN_MODELS_DIR:-$DATA_REPO/Models/gapseq/models}"
TARA_MEDIUM="$DATA_REPO/Models/gapseq/media/mathomics.txt"
BUNDLED_MEDIUM="$REPO_ROOT/experiments/ocean_samples/media/ocean_complete_media.csv"
DEFAULT_MATRIX="$DATA_REPO/Data/CEODOS_MAG_TPMs_IDs_matrix.tsv"
MATRIX_INPUT="${1:-${OCEAN_ABUNDANCE_MATRIX:-$DEFAULT_MATRIX}}"

PRESENCE_THRESHOLD="${OCEAN_PRESENCE_THRESHOLD:-0}"
ALPHA="${OCEAN_ALPHA:-0.20}"
MINIMAL_GROWTH="${OCEAN_MINIMAL_GROWTH:-0.01}"
PRODUCT_RETENTION="${OCEAN_PRODUCT_RETENTION:-0.90}"
LEX_TOL="${OCEAN_LEX_TOL:-1e-5}"
PROD_TOL="${OCEAN_PROD_TOL:-1e-6}"
INT_TOL="${OCEAN_INT_TOL:-1e-9}"
UPTAKE_BOUND="${OCEAN_UPTAKE_BOUND:--1000}"
MAX_MIN_COMMS="${OCEAN_MAX_MIN_COMMS:-100}"
MAX_MAGS_PER_SAMPLE="${OCEAN_MAX_MAGS_PER_SAMPLE:-}"
ALLOW_MISSING_MODELS="${OCEAN_ALLOW_MISSING_MODELS:-0}"

PARTITION="${OCEAN_PARTITION:-general}"
TIME_LIMIT="${OCEAN_TIME_LIMIT:-12:00:00}"
CPUS_PER_TASK="${OCEAN_CPUS_PER_TASK:-4}"
MEMORY="${OCEAN_MEMORY:-32G}"
MAX_CONCURRENT="${OCEAN_MAX_CONCURRENT:-2}"
WLS_SESSIONS="${OCEAN_WLS_SESSIONS:-2}"

MODEL_MAP="${OCEAN_MODEL_MAP:-}"

if (( MAX_CONCURRENT > WLS_SESSIONS )); then
    echo "ERROR: OCEAN_MAX_CONCURRENT=$MAX_CONCURRENT exceeds the configured" >&2
    echo "Gurobi WLS session limit OCEAN_WLS_SESSIONS=$WLS_SESSIONS." >&2
    echo "Increase OCEAN_WLS_SESSIONS only if the license actually supports it." >&2
    exit 1
fi

if [[ ! -d "$DATA_REPO" ]]; then
    echo "Data repository not found. Cloning to: $DATA_REPO"
    git clone --depth 1 \
        https://github.com/mathomics/tara_chile_metabolic_models.git \
        "$DATA_REPO"
fi

if [[ -n "${OCEAN_MEDIUM_FILE:-}" ]]; then
    MEDIUM_FILE="$OCEAN_MEDIUM_FILE"
elif [[ -f "$TARA_MEDIUM" ]]; then
    MEDIUM_FILE="$TARA_MEDIUM"
else
    MEDIUM_FILE="$BUNDLED_MEDIUM"
    echo "mathomics.txt not found; using bundled temporary medium:"
    echo "  $MEDIUM_FILE"
fi

MATRIX_PATH="$(readlink -f "$MATRIX_INPUT")"

if [[ ! -f "$MATRIX_PATH" ]]; then
    echo "ERROR: abundance matrix not found: $MATRIX_PATH" >&2
    echo "Expected by default:" >&2
    echo "  $DATA_REPO/Data/CEODOS_MAG_TPMs_IDs_matrix.tsv" >&2
    exit 1
fi

if [[ ! -d "$MODELS_DIR" ]]; then
    echo "ERROR: gapseq models directory not found: $MODELS_DIR" >&2
    exit 1
fi

if [[ ! -f "$MEDIUM_FILE" ]]; then
    echo "ERROR: medium file not found: $MEDIUM_FILE" >&2
    echo "Set OCEAN_MEDIUM_FILE explicitly or restore the bundled medium." >&2
    exit 1
fi

RUN_STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="${OCEAN_RUN_ROOT:-$HOME/misosoup_runs/ocean_samples_${RUN_STAMP}}"
MANIFEST="$RUN_ROOT/manifest.yaml"

mkdir -p "$RUN_ROOT/logs" "$RUN_ROOT/samples" "$RUN_ROOT/work"

echo "============================================================"
echo "Preparing ocean-sample run"
echo "MiSoSoup repo       : $REPO_ROOT"
echo "Abundance matrix    : $MATRIX_PATH"
echo "TARA repo           : $DATA_REPO"
echo "Models              : $MODELS_DIR"
echo "Medium              : $MEDIUM_FILE"
echo "Run root            : $RUN_ROOT"
echo "Presence threshold  : $PRESENCE_THRESHOLD"
echo "Max MAGs/sample     : ${MAX_MAGS_PER_SAMPLE:-all}"
echo "Allow missing models: $ALLOW_MISSING_MODELS"
echo "Partition           : $PARTITION"
echo "Time limit/task     : $TIME_LIMIT"
echo "CPUs/task           : $CPUS_PER_TASK"
echo "Memory/task         : $MEMORY"
echo "Array concurrency   : $MAX_CONCURRENT"
echo "============================================================"

PREPARE_ARGS=(
    --matrix "$MATRIX_PATH"
    --models-dir "$MODELS_DIR"
    --output "$MANIFEST"
    --presence-threshold "$PRESENCE_THRESHOLD"
    --expected-samples 159
    --expected-mags 1375
)

if [[ -n "$MAX_MAGS_PER_SAMPLE" ]]; then
    PREPARE_ARGS+=(--max-mags-per-sample "$MAX_MAGS_PER_SAMPLE")
fi

if [[ "$ALLOW_MISSING_MODELS" == "1" ]]; then
    PREPARE_ARGS+=(--allow-missing-models)
fi

if [[ -n "$MODEL_MAP" ]]; then
    PREPARE_ARGS+=(--model-map "$MODEL_MAP")
fi

python -u "$SCRIPT_DIR/prepare_manifest.py" "${PREPARE_ARGS[@]}"

N_SAMPLES="$(
    python - "$MANIFEST" <<'PY'
import sys
import yaml
with open(sys.argv[1], "r", encoding="utf-8") as handle:
    manifest = yaml.safe_load(handle)
print(len(manifest["samples"]))
PY
)"

if [[ "$N_SAMPLES" -lt 1 ]]; then
    echo "ERROR: manifest contains no samples." >&2
    exit 1
fi

ARRAY_LAST=$((N_SAMPLES - 1))
ARRAY_RANGE="${OCEAN_ARRAY_RANGE:-0-${ARRAY_LAST}}"
ARRAY_SPEC="${ARRAY_RANGE}%${MAX_CONCURRENT}"

export OCEAN_MANIFEST="$MANIFEST"
export OCEAN_OUTPUT_ROOT="$RUN_ROOT"
export OCEAN_MEDIUM_FILE="$MEDIUM_FILE"
export OCEAN_ALPHA="$ALPHA"
export OCEAN_MINIMAL_GROWTH="$MINIMAL_GROWTH"
export OCEAN_PRODUCT_RETENTION="$PRODUCT_RETENTION"
export OCEAN_LEX_TOL="$LEX_TOL"
export OCEAN_PROD_TOL="$PROD_TOL"
export OCEAN_INT_TOL="$INT_TOL"
export OCEAN_UPTAKE_BOUND="$UPTAKE_BOUND"
export OCEAN_MAX_MIN_COMMS="$MAX_MIN_COMMS"

ARRAY_SUBMIT="$(
    sbatch --parsable \
        --partition="$PARTITION" \
        --time="$TIME_LIMIT" \
        --cpus-per-task="$CPUS_PER_TASK" \
        --mem="$MEMORY" \
        --array="$ARRAY_SPEC" \
        --output="$RUN_ROOT/logs/ocean_%A_%a.out" \
        --error="$RUN_ROOT/logs/ocean_%A_%a.err" \
        "$SCRIPT_DIR/sample_array.slurm"
)"
ARRAY_JOB_ID="${ARRAY_SUBMIT%%;*}"

COLLECT_SUBMIT="$(
    sbatch --parsable \
        --dependency="afterany:${ARRAY_JOB_ID}" \
        --partition=debug \
        --time=00:10:00 \
        --cpus-per-task=1 \
        --mem=2G \
        --output="$RUN_ROOT/logs/collect_%j.out" \
        --error="$RUN_ROOT/logs/collect_%j.err" \
        "$SCRIPT_DIR/collect_results.slurm"
)"
COLLECT_JOB_ID="${COLLECT_SUBMIT%%;*}"

cat > "$RUN_ROOT/submission_info.txt" <<EOF
run_root=$RUN_ROOT
manifest=$MANIFEST
matrix=$MATRIX_PATH
models_dir=$MODELS_DIR
medium_file=$MEDIUM_FILE
number_samples=$N_SAMPLES
array_job_id=$ARRAY_JOB_ID
collector_job_id=$COLLECT_JOB_ID
array_spec=$ARRAY_SPEC
partition=$PARTITION
time_limit=$TIME_LIMIT
cpus_per_task=$CPUS_PER_TASK
memory=$MEMORY
max_concurrent=$MAX_CONCURRENT
max_mags_per_sample=$MAX_MAGS_PER_SAMPLE
allow_missing_models=$ALLOW_MISSING_MODELS
alpha=$ALPHA
minimal_growth=$MINIMAL_GROWTH
product_retention=$PRODUCT_RETENTION
presence_threshold=$PRESENCE_THRESHOLD
uptake_bound=$UPTAKE_BOUND
EOF

echo
echo "============================================================"
echo "Submitted ocean-sample workflow"
echo "============================================================"
echo "Array job      : $ARRAY_JOB_ID"
echo "Collector job  : $COLLECT_JOB_ID"
echo "Array          : $ARRAY_SPEC"
echo "Run root       : $RUN_ROOT"
echo
echo "Status:"
echo "  squeue -j $ARRAY_JOB_ID,$COLLECT_JOB_ID"
echo
echo "Accounting:"
echo "  sacct -j $ARRAY_JOB_ID --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode"
echo
echo "Progress:"
echo "  find \"$RUN_ROOT/samples\" -name result.yaml | wc -l"
echo "  find \"$RUN_ROOT/samples\" -name failure.yaml | wc -l"
echo
echo "After collector finishes:"
echo "  column -s, -t < \"$RUN_ROOT/sample_summary.csv\" | less -S"
echo "============================================================"
