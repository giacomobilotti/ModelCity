# Generated figures

Outputs are organized first by dataset and then by analysis:

```text
figures/
├── yautepec/
│   ├── unified/            system metrics and plots
│   ├── sampling-python/    chronological sampling
│   ├── longevity-python/   longevity
│   ├── sampling/           R sampling reference
│   ├── longevity/          R longevity reference
│   └── analysis/           dataset-specific R analysis
└── viabundus/
    ├── unified/            system metrics and plots
    ├── sampling-python/    chronological sampling
    ├── longevity-python/   longevity
    └── analysis/           dataset-specific R analysis
```

Generate the three Python analyses with:

```bash
python3 scripts/run_unified.py
python3 scripts/run_sampling.py
python3 scripts/run_longevity.py
python3 scripts/run_longevity.py --dataset viabundus
```

`run_unified.py` and `run_sampling.py` run every example dataset unless
`--dataset` is set. `run_longevity.py` runs Yautepec unless `--dataset` names
another survey.

Generate the dataset-specific `analysis/` directories with:

```bash
Rscript scripts/r/run-plots.R
Rscript scripts/r/run-viabundus-plots.R
```

The Yautepec `sampling/` and `longevity/` directories are the R reference
outputs, from `scripts/r/run-sampling-plots.R` and
`scripts/r/run-longevity-plots.R`.

Figures and tables are generated artifacts. Do not place dataset outputs
directly in this directory; add them beneath
`figures/<dataset>/<workflow>/`.
