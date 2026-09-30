
# Ocean-sample product and minimum-community workflow

This experiment scales the validated MiSoSoup product-preservation workflow to
the 159 Chilean ocean metagenomic samples.

Each sample is defined by the MAGs with positive abundance in a 159 x 1375
sample-by-MAG abundance matrix. The abundance value is used **only to decide
whether a MAG is present**; it is not currently used as a flux weight or biomass
fraction.

For every sample, the workflow:

1. loads only the MAG models present in that sample;
2. builds the full sample community;
3. applies the headerless `mathomics.txt` medium;
4. computes individual maximum secretion for every producible global exchange;
5. applies the C/N/S/P product filter;
6. runs Stage A to obtain `K*`, the maximum number of simultaneous products;
7. runs Stage B to obtain `Q*` and the selected product identities;
8. runs the validated fresh Stage C with finite split-flux bounds and pFBA;
9. runs Stage D to obtain the globally minimum community size and enumerate
   alternative minimum communities.

The 159 samples are submitted as a Slurm **job array**. A second dependent job
collects all results into analysis-ready CSV files.

## Input data

Default data repository:

```text
~/tara_chile_metabolic_models
```

Expected model directory:

```text
Models/gapseq/models/
```

Final expected medium:

```text
Models/gapseq/media/mathomics.txt
```

Until that file is added to the TARA Chile repository, the workflow falls back
automatically to the temporary medium bundled in MiSoSoup:

```text
experiments/ocean_samples/media/complex_media_gapseq2.csv
```

The temporary CSV contains 125 compounds with columns `compound`, `name`, and
`maxFlux`. Its `maxFlux` values are used directly as uptake magnitudes; for
example `maxFlux = 10` becomes a global-exchange lower bound of `-10`.

The future `mathomics.txt` format is also supported: one metabolite per row and
no header. Headerless entries use the fallback uptake lower bound `-1000`, which
can be changed with `OCEAN_UPTAKE_BOUND`.

The abundance matrix is already available in the TARA Chile repository at:

```text
Data/CEODOS_MAG_TPMs_IDs_matrix.tsv
```

The launcher now uses that file by default. It must contain one header row with
1375 MAG identifiers and one leading column containing the 159 sample
identifiers. CSV and TSV are supported. If the matrix is accidentally
transposed, the preparatory step will detect a 1375 x 159 matrix and transpose
it automatically.

## MAG-to-model matching

The preparatory step scans all `.xml` and `.sbml` files below the gapseq model
directory and matches abundance-matrix MAG identifiers to model filenames.

Matching is conservative:

1. exact filename/stem match;
2. normalized match after removing common `.xml`, `.sbml`, `model`, and
   `gapseq` decorations.

Every present MAG must map to exactly one model and one model may not silently
represent two different MAG columns.

If the names do not match automatically, provide an explicit two-column table:

```text
mag_id    model_path
MAG_001   /path/to/MAG_001.xml
MAG_002   /path/to/MAG_002.xml
```

and export:

```bash
export OCEAN_MODEL_MAP=/path/to/mag_model_map.tsv
```

before submission.

## Medium matching

For each sample, medium tokens are resolved against the sample's global
exchange reactions using several common representations, including:

```text
R_EX_glc__D_e
EX_glc__D_e
M_glc__D_e
glc__D_e
glc__D
```

and gapseq-style external suffixes such as `_e0`.

Medium metabolites that are absent from a particular sample are recorded as
unresolved for that sample but do not cause an error. Ambiguous matches do
cause an error.

## Parallelization and the Gurobi license

The current Gurobi WLS academic license allows two concurrent sessions.
Therefore the launcher uses:

```text
0-158%2
```

by default: 159 array tasks, with at most two sample analyses running at once.

Do **not** raise `OCEAN_MAX_CONCURRENT` above 2 unless the Gurobi license limit
has changed. Otherwise tasks can fail before optimization starts.

The default production partition is `general`, not `debug`. NLHPC documents
`debug` as a short test partition, whereas the full 159-sample workflow may
need longer jobs.

Default resources per sample:

```text
partition = general
CPUs      = 4
memory    = 32G
time      = 12:00:00
```

These can be overridden at submission time.

Whether the whole analysis finishes within 12 hours depends on the MAG count
and optimization time of the real samples. With two concurrent Gurobi
sessions, a 12-hour wall-clock target requires an average sample runtime of
roughly 9 minutes or less. The output records per-sample timings so the
projection can be checked after the first few tasks finish.

## Submit the complete workflow

First update MiSoSoup:

```bash
cd ~/misosoup_product
git pull origin main
```

The repository defaults are now sufficient even before `mathomics.txt` is uploaded:

```bash
bash experiments/ocean_samples/submit_ocean_samples.sh
```

This uses:

```text
~/tara_chile_metabolic_models/Data/CEODOS_MAG_TPMs_IDs_matrix.tsv
~/tara_chile_metabolic_models/Models/gapseq/models/
~/tara_chile_metabolic_models/Models/gapseq/media/mathomics.txt  # preferred when present
```

You can still override the matrix or repository path explicitly if needed.

If the data repository does not yet exist at that path, the launcher attempts a
shallow clone from:

```text
https://github.com/mathomics/tara_chile_metabolic_models.git
```

Existing clones are never automatically modified.

The launcher validates the matrix and MAG/model mapping, creates a run
manifest, submits the 159-element array, limits it to two simultaneous Gurobi
processes, and submits a collector with `afterany` dependency so aggregation
runs even if individual samples fail.

## Recommended pilot before the 159-sample run

Because the abundance matrix and `mathomics.txt` were not available while this
workflow was implemented, first run a small pilot after those files appear:

```bash
export OCEAN_ARRAY_RANGE=0-3
bash experiments/ocean_samples/submit_ocean_samples.sh
```

The concurrency cap is still applied, so this becomes `0-3%2`. The collector
will mark the other samples as missing in this pilot run; that is expected.

After confirming model-name matching, medium matching, feasibility, memory, and
per-sample runtime, unset the pilot range and submit all 159 samples:

```bash
unset OCEAN_ARRAY_RANGE
```

## Useful overrides

```bash
export OCEAN_PARTITION=main
export OCEAN_TIME_LIMIT=06:00:00
export OCEAN_CPUS_PER_TASK=4
export OCEAN_MEMORY=48G
export OCEAN_MAX_CONCURRENT=2
export OCEAN_PRESENCE_THRESHOLD=0
export OCEAN_ALPHA=0.20
export OCEAN_MINIMAL_GROWTH=0.01
export OCEAN_PRODUCT_RETENTION=0.90
export OCEAN_UPTAKE_BOUND=-1000
export OCEAN_MAX_MIN_COMMS=100
```

A MAG is considered present when `abundance > OCEAN_PRESENCE_THRESHOLD`; the
default threshold is `0`.

## Per-sample output

A run is created below:

```text
~/misosoup_runs/ocean_samples_YYYYMMDD_HHMMSS/
```

Per-sample results are stored as `samples/000/result.yaml` through
`samples/158/result.yaml`. A failed task instead creates `failure.yaml` and
`traceback.txt` in its sample directory.

Each successful result contains the sample MAGs and abundances, medium audit,
individual product maxima, Stage-A `K*`, Stage-B `Q*`, Stage-C reference
products/fluxes and pFBA objective, Stage-D minimum community size and
communities, and per-stage timing.

## Aggregated output

After all array elements terminate, the collector writes:

```text
sample_summary.csv
producible_products_long.csv
selected_products_long.csv
minimal_communities.csv
minimal_community_members_long.csv
aggregate.yaml
```

`producible_products_long.csv` contains all C/N/S/P-filtered individually
producible products and their sample-specific maximum secretion.
`selected_products_long.csv` contains the final Stage-C reference products and
fluxes. `minimal_communities.csv` contains one row per minimum community, and
`minimal_community_members_long.csv` contains one row per MAG membership.

## Monitor the array

```bash
squeue -j <ARRAY_JOB_ID>,<COLLECTOR_JOB_ID>
```

```bash
sacct -j <ARRAY_JOB_ID> \
  --format=JobID,JobName,State,Elapsed,MaxRSS,ExitCode
```

Progress:

```bash
find <RUN_ROOT>/samples -name result.yaml | wc -l
find <RUN_ROOT>/samples -name failure.yaml | wc -l
```

## Important interpretation

The abundance matrix currently determines sample membership only. A MAG with
abundance `0.0001` and one with abundance `0.5` are both represented once if
both are above the presence threshold. An abundance-weighted community
formulation would be a distinct model and should be evaluated separately rather
than introduced silently.
