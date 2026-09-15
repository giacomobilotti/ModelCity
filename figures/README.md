# Generated figures

Outputs are organized first by dataset and then by workflow:

```text
figures/
├── yautepec/
│   ├── unified/    shared ModelCity plots and metric tables
│   └── analysis/   dataset-specific R analysis outputs
└── viabundus/
    ├── unified/    shared ModelCity plots and metric tables
    └── analysis/   dataset-specific R analysis outputs
```

Generate the `unified/` directories with:

```bash
python3 scripts/run_unified.py
```

Generate the dataset-specific `analysis/` directories with:

```bash
Rscript scripts/r/run-plots.R
Rscript scripts/r/run-viabundus-plots.R
```

Figures and tables are generated artifacts. Do not place dataset outputs
directly in this directory; add them beneath
`figures/<dataset>/<workflow>/`.
