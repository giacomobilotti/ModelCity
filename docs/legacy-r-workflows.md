# Legacy R and Quarto workflows

The reusable ModelCity interface is the Python package under `src/modelcity/`.
The R runners under [`scripts/r/`](../scripts/r/) and Quarto files under
[`notebooks/`](../notebooks/) are dataset-specific project analyses. They are
retained for provenance and for figures that are outside the shared Python
metric set.

They should not be used as templates for adapting a new dataset.

## Scope

| File | Purpose | Status |
| --- | --- | --- |
| `scripts/r/run-plots.R` | CSV-based Yautepec plots, including trajectories and Sankey figures | Standalone replacement for parts of the original notebook |
| `scripts/r/run-viabundus-plots.R` | Viabundus maps, regional summaries, alluvial plots, and rank trajectories | Uses approximations where upstream assets are absent |
| `notebooks/master-yautepec.qmd` | Original Yautepec survey manuscript analysis | Requires spatial data and assets not included here |
| `notebooks/analyses.qmd` | Original Viabundus paper analysis | Requires upstream scripts and spatial, climate, and linguistic data |

## CSV-based R runners

`run-plots.R` reads `data/examples/yautepec.csv` and writes dataset-specific
outputs to `figures/yautepec/analysis/`. It uses packages including:

- `tidyverse`;
- `ggplot2`;
- `ggsankey`;
- `rcarbon`.

`run-viabundus-plots.R` reads `data/examples/viabundus.csv` and writes to
`figures/viabundus/analysis/`. It uses packages including:

- `tidyverse`;
- `sf`;
- `ggalluvial`;
- `patchwork`;
- `scales`;
- optionally `rnaturalearth` for map backgrounds.

Both runners derive the repository root from their own location when invoked
with `Rscript`, and use a project-local `R_libs/` directory when it exists:

```bash
Rscript scripts/r/run-plots.R
Rscript scripts/r/run-viabundus-plots.R
```

No `renv` lockfile or complete automated R environment setup is provided.

The Viabundus runner replaces the unavailable
`linguistic_areas.gpkg` with approximate coordinate rules. Those regions are
suitable for reproducing the included grouping figures, not for substantive
claims about historical language boundaries.

## Original Quarto notebooks

`notebooks/master-yautepec.qmd` expects:

- `data/derived_data/yautepec.gpkg` with several layers;
- `survey_sites.png`;
- `references.bib`;
- a broader R geospatial and statistical environment.

`notebooks/analyses.qmd` expects resources including:

- `analyses/00_loading_data.R`;
- `data/raw_data/NOAA`;
- `data/derived_data/linguistic_areas.gpkg`;
- map data and several R spatial/plotting packages.

These dependencies are not all present in the repository, so a clean checkout
cannot render the notebooks in full. The Markdown documentation must not imply
otherwise.

## Relationship to the Python toolkit

Some R calculations motivated or validate shared Python behavior, including:

- absolute KDE bandwidth compatible with R's `density(..., bw=...)`;
- BP and BC/AD conversion parity;
- settlement duration, event, turnover, size, rank-size, and threshold
  summaries.

The Python demonstration runner
[`scripts/run_unified.py`](../scripts/run_unified.py) is the maintained example
of applying one common analysis to different schemas. R-only maps, alluvial
diagrams, Sankey figures, and paper-specific models are not currently part of
the general Python API.

For a new dataset, begin with [Adapting a new dataset](adapting-datasets.md),
not these scripts.
