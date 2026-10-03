#!/bin/bash
set -e
set -u
set -o pipefail

REPO_ROOT="${JACCARD_REPO_ROOT:-$HOME/misosoup_product}"
DATA_ROOT="${JACCARD_DATA_ROOT:-$HOME/tara_chile_metabolic_models}"
PARTITION="${JACCARD_PARTITION:-general}"
TIME_LIMIT="${JACCARD_TIME_LIMIT:-03:00:00}"
MEMORY="${JACCARD_MEMORY:-16G}"

STAMP="$(date +%Y%m%d_%H%M%S)"
RUN_ROOT="${JACCARD_OUTPUT_ROOT:-$HOME/misosoup_runs/tara_jaccard_${STAMP}}"
mkdir -p "$RUN_ROOT/logs"

echo "Validating the three 1,375-model sets before submission..."
python -u "$REPO_ROOT/experiments/ocean_samples/jaccard_models.py"   --data-root "$DATA_ROOT"   --output-root "$RUN_ROOT"   --expected-models 1375   --validate-only

JOB_ID="$(
  cd "$RUN_ROOT/logs"
  sbatch --parsable     --partition="$PARTITION"     --time="$TIME_LIMIT"     --mem="$MEMORY"     --export=ALL,JACCARD_DATA_ROOT="$DATA_ROOT",JACCARD_OUTPUT_ROOT="$RUN_ROOT"     "$REPO_ROOT/experiments/ocean_samples/jaccard_models.slurm"
)"

cat <<EOF
Submitted TARA Chile reaction-Jaccard job.

Job ID   : $JOB_ID
Run root : $RUN_ROOT

Monitor:
  squeue -j $JOB_ID

Accounting:
  sacct -j $JOB_ID --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode

Live log:
  tail -f $RUN_ROOT/logs/tara_jaccard_${JOB_ID}.out

Outputs (one directory per model set):
  $RUN_ROOT/draft/
  $RUN_ROOT/complete_medium/
  $RUN_ROOT/LS2N_medium/

Each contains:
  jaccard_reaction_distance.npy
  jaccard_reaction_distance_upper.npy
  jaccard_reaction_distance.npz
  jaccard_reaction_distance.tsv.gz
  models.tsv
  summary.json
EOF
