# ModelCity

ModelCity is a Python toolkit for comparing how settlement systems change over
long periods of time. It maps differently structured tables onto one canonical
data model, then applies the same three analyses whether settlement size is
recorded as occupied area or population: system metrics and plots, chronological
Monte Carlo sampling, and longevity.

The package is dataset-agnostic. The Yautepec and Viabundus files in this
repository are worked examples, not requirements.

ModelCity was developed within the ERC project *The Model City: Drivers and
Mechanisms of Long-term Urban Evolution and Resilience* (PI: Iza Romanowska).

## How it works

```text
your DataFrame
    → map id, size, and time columns
    → canonical Settlements object
         ├── system metrics and plots
         ├── chronological sampling
         └── longevity
    → pandas tables and Matplotlib figures
```

The three analyses share that object:

- **System metrics and plots** summarise duration, foundations and
  abandonments, churn, size, rank-size, and large-centre share.
- **Chronological sampling** draws foundation, peak, and abandonment dates
  inside phase bounds and carries each draw through growth and persistence.
- **Longevity** treats a settlement still occupied at the end of the record as
  right-censored, fits a survival model to the sampled durations, and imputes
  a total lifespan.

ModelCity supports:

- interval records with explicit start and end years, when settlements use a
  shared set of period bounds;
- snapshot records with one observation year;
- sizes measured as area or population;
- source dates in years BP or calendar years;
- BP or BC/AD labels without changing the underlying calculations.

Current plot labels use “site” for area data and “city” for population data.
These are presentation conventions, not restrictions on which settlements can
be analyzed.

## Installation

ModelCity requires Python 3.9 or newer. From the repository root:

```bash
python3 -m pip install -e .
```

Install the test tools with:

```bash
python3 -m pip install -e ".[dev]"
```

The optional `geo` extra installs `pyproj`, which the Viabundus loader in
`scripts/run_unified.py` uses for region assignment. The sampling and
longevity runners import that loader:

```bash
python3 -m pip install -e ".[dev,geo]"
```

## Quick start with your own data

The example below uses an interval dataset whose dates are in years BP:

```python
import pandas as pd
import modelcity as mc

frame = pd.DataFrame(
    {
        "place": ["A", "A", "B", "B"],
        "extent_ha": [12.0, 18.0, 7.0, 0.0],
        "from_bp": [3200, 2900, 3200, 2900],
        "to_bp": [2900, 2600, 2900, 2600],
    }
)

settlements = mc.as_settlements(
    frame,
    id="place",
    size="extent_ha",
    start="from_bp",
    end="to_bp",
    size_type="area",
    time_scale="BP",
)

totals = mc.size_totals(settlements)
figure = mc.plot_size_totals(settlements)
mc.save_figure(figure, "outputs/size-totals.png")
```

For snapshot observations, omit `end`. ModelCity derives contiguous period
bounds from the midpoints between observation years:

```python
snapshots_frame = pd.DataFrame(
    {
        "city_id": ["A", "A", "B", "B"],
        "population": [900, 1200, 500, 700],
        "year": [1200, 1300, 1200, 1300],
    }
)

snapshots = mc.as_settlements(
    snapshots_frame,
    id="city_id",
    size="population",
    start="year",
    size_type="population",
    time_scale="CE",
)
```

See [Adapting datasets](docs/adapting-datasets.md) for complete mapping recipes.

## Documentation

- [Documentation home](docs/index.md)
- [Getting started](docs/getting-started.md)
- [Canonical data model](docs/data-model.md)
- [End-to-end workflow](docs/workflow.md)
- [Methodologies by stage](docs/methodologies.md)
- [Metrics and plots](docs/metrics-and-plots.md) (`scripts/run_unified.py`)
- [Chronological sampling](docs/chronological-sampling.md) (`scripts/run_sampling.py`)
- [Longevity](docs/longevity.md) (`scripts/run_longevity.py`)
- [Adapting new datasets](docs/adapting-datasets.md)
- [Python API](docs/python-api.md)
- [Reproducibility](docs/reproducibility.md)
- [Legacy R and Quarto workflows](docs/legacy-r-workflows.md)

## Included demonstrations

Source-specific loaders for the example CSVs live in
[`scripts/run_unified.py`](scripts/run_unified.py), outside the package. The
sampling and longevity runners import those loaders. Each runner then applies
one of the three analyses:

| Runner | Analysis | Output |
| --- | --- | --- |
| [`scripts/run_unified.py`](scripts/run_unified.py) | System metrics and plots | `figures/<dataset>/unified/` |
| [`scripts/run_sampling.py`](scripts/run_sampling.py) | Chronological sampling | `figures/<dataset>/sampling-python/` |
| [`scripts/run_longevity.py`](scripts/run_longevity.py) | Longevity | `figures/<dataset>/longevity-python/` |

```bash
python3 scripts/run_unified.py
python3 scripts/run_sampling.py
python3 scripts/run_longevity.py
python3 scripts/run_longevity.py --dataset viabundus
```

`run_unified.py` and `run_sampling.py` run every example dataset unless
`--dataset` is set. `run_longevity.py` runs Yautepec unless `--dataset` names
another survey. Pass `--dataset` to any of them to select one survey.

## Repository layout

```text
src/modelcity/        reusable Python package
tests/                unit and integration tests
docs/                 user and developer guides
data/examples/         included example datasets
scripts/               runners for the three analyses, plus dataset-specific R runners
notebooks/             project-specific Quarto analyses
figures/<dataset>/     unified, sampling, longevity, and project-analysis outputs
```

## Testing

```bash
pytest
```

The tests cover schema validation, interval and snapshot normalization, time
conversion, metrics, plotting, chronological sampling, longevity, and
integration with the included examples.

## Citation, license, and contact

TBA
