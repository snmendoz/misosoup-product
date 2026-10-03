# Individual growth QC for the ocean MAG models

This workflow checks each MAG metabolic model independently before it is used
inside a microbial community.

For every model it performs two growth maximizations:

1. **Selected medium**: the current ocean medium, currently the bundled
   `complex_media_gapseq2.csv` until `mathomics.txt` is available.
2. **Rich medium**: every global exchange reaction in the one-member community
   is opened for uptake with lower bound `-1000` by default.

The model is classified as:

- `grows_in_medium`: maximum growth reaches the requested threshold in the
  selected medium.
- `medium_limited`: it fails in the selected medium but reaches the threshold
  when all global exchanges are open.
- `no_growth_even_rich`: even the rich-medium optimization remains below the
  threshold.
- `solver_error`: the model could not be evaluated successfully.

The default growth threshold is `0.01`.

## Immediate diagnostic for MMC_B0017 and MMC_B0965

```bash
cd ~/misosoup_product
git pull origin main

export OCEAN_QC_MAG_IDS=MMC_B0017,MMC_B0965

bash experiments/ocean_samples/model_qc/submit_model_qc.sh
```

This creates a two-model manifest and normally submits a single QC chunk.

## Full 1375-model QC

```bash
unset OCEAN_QC_MAG_IDS

bash experiments/ocean_samples/model_qc/submit_model_qc.sh
```

By default the 1375 models are divided into chunks of 25 models. A maximum of
two chunks run simultaneously because the current Gurobi WLS setup supports two
concurrent sessions.

The final outputs are:

```text
model_qc.csv
aggregate.yaml
models/0000/result.yaml
models/0001/result.yaml
...
```

The CSV contains the MAG ID, biomass reaction, maximum growth in the selected
medium, maximum growth with all exchanges open, number of selected-medium
metabolites matched, and the diagnostic class.
