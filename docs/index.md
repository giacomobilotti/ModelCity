# ModelCity documentation

ModelCity is a dataset-agnostic Python toolkit for normalizing longitudinal
settlement tables and applying three shared analyses: system metrics and
plots, chronological Monte Carlo sampling, and longevity.

## Start here

- New user: [Getting started](getting-started.md)
- Bringing a table: [Adapting a new dataset](adapting-datasets.md)
- Understanding assumptions: [Canonical data model](data-model.md)
- Building an analysis: [End-to-end workflow](workflow.md)
- Library stages and functions: [Methodologies by stage](methodologies.md)

## The three analyses

Each analysis has a guide and a runner. The runners share the example loaders
and write under `figures/<dataset>/`.

- [Metrics and plots](metrics-and-plots.md) —
  [`scripts/run_unified.py`](../scripts/run_unified.py) writes
  `figures/<dataset>/unified/`.
- [Chronological sampling](chronological-sampling.md) —
  [`scripts/run_sampling.py`](../scripts/run_sampling.py) writes
  `figures/<dataset>/sampling-python/`.
- [Longevity](longevity.md) —
  [`scripts/run_longevity.py`](../scripts/run_longevity.py) writes
  `figures/<dataset>/longevity-python/`.

## Reference

- [Python API](python-api.md)
- [Reproducibility](reproducibility.md)

## Project-specific material

- [Legacy R and Quarto workflows](legacy-r-workflows.md)

The Yautepec and Viabundus datasets demonstrate the adapter pattern used by
all three runners. Any table that maps onto the canonical model can use the
same analyses. The package runs on a table you supply.
