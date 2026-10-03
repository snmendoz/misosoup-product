# Pairwise reaction-content Jaccard distances

This analysis compares all 1,375 TARA Chile gapseq models within each of three
model sets:

- `Models/gapseq/models/draft`
- `Models/gapseq/models/gapfilled/complete_medium`
- `Models/gapseq/models/gapfilled/LS2N_medium`

For models (i) and (j), let (R_i) and (R_j) be their sets of exact SBML
reaction IDs. The reported distance is

[
d_J(i,j) = 1 - \frac{|R_i \cap R_j|}{|R_i \cup R_j|}.
]

The implementation evaluates only (i<j), i.e. 944,625 unique pairs for 1,375
models. It then mirrors the values to return the requested full symmetric
1,375 x 1,375 matrix with a zero diagonal.

Before computing anything, the script requires exactly 1,375 SBML files in
each set and verifies that the normalized model IDs are identical across all
three sets. Draft filenames such as `MMC_A0001-draft.xml` are normalized to
`MMC_A0001`, so all three matrices use exactly the same row/column order.

## Run on Leftraru

```bash
cd ~/misosoup_product
git pull origin main
conda activate misosoup

bash experiments/ocean_samples/submit_jaccard_models.sh
```

The launcher validates the inputs before calling `sbatch`.

Default output:

```text
~/misosoup_runs/tara_jaccard_YYYYMMDD_HHMMSS/
├── model_ids.tsv
├── run_summary.json
├── draft/
├── complete_medium/
└── LS2N_medium/
```

Each model-set directory contains:

- `jaccard_reaction_distance.npy`: full float32 symmetric matrix.
- `jaccard_reaction_distance_upper.npy`: the 944,625 computed (i<j)
  distances in row-major upper-triangle order.
- `jaccard_reaction_distance.npz`: compressed matrix plus model IDs.
- `jaccard_reaction_distance.tsv.gz`: labeled, compressed text matrix.
- `models.tsv`: matrix index, normalized model ID, reaction count and source
  SBML path.
- `summary.json`: reaction-count and Jaccard-distance summaries.

The calculation does not use FBA, Gurobi or any medium constraints. It compares
reaction content only.

## Useful overrides

```bash
export JACCARD_PARTITION=general
export JACCARD_TIME_LIMIT=03:00:00
export JACCARD_MEMORY=16G
export JACCARD_DATA_ROOT=~/tara_chile_metabolic_models
bash experiments/ocean_samples/submit_jaccard_models.sh
```
