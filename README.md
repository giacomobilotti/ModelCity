# ModelCity

ModelCity is a Python toolkit for comparing how settlement systems change over
long periods of time. It maps differently structured tables onto one canonical
data model, then applies the same metrics and plots whether settlement size is
recorded as occupied area or population.

The package is dataset-agnostic. The Yautepec and Viabundus files in this
repository are worked examples, not requirements.

ModelCity was developed within the ERC project *The Model City: Drivers and
Mechanisms of Long-term Urban Evolution and Resilience* (PI: Iza Romanowska).

## How it works

```text
your DataFrame
    → map id, size, and time columns
    → canonical Settlements object
    → durations, events, churn, size summaries, rank-size, and urban share
    → pandas tables and Matplotlib figures
```

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

The optional `geo` extra installs `pyproj`, which the Viabundus demo loader in
`scripts/run_unified.py` uses for region assignment:

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
- [Metrics and plots](docs/metrics-and-plots.md)
- [Adapting new datasets](docs/adapting-datasets.md)
- [Python API](docs/python-api.md)
- [Reproducibility](docs/reproducibility.md)
- [Legacy R and Quarto workflows](docs/legacy-r-workflows.md)

## Included demonstrations

Source-specific loaders for the example CSVs live in
[`scripts/run_unified.py`](scripts/run_unified.py), not in the package. They
map each table through `as_settlements()`, then the script applies the same
seven plotting calls to both datasets.

```bash
python3 scripts/run_unified.py
python3 scripts/run_unified.py --dataset yautepec
```

Figures and tables are written under `figures/<dataset>/unified/`.

## Repository layout

```text
src/modelcity/        reusable Python package
tests/                unit and integration tests
docs/                 user and developer guides
data/examples/         included example datasets
scripts/               Python runner and dataset-specific R runners
notebooks/             project-specific Quarto analyses
figures/<dataset>/     unified and project-analysis outputs
```

## Testing

```bash
pytest
```

The tests cover schema validation, interval and snapshot normalization, time
conversion, metrics, plotting, and integration with the included examples.

## Citation, license, and contact

TBA
