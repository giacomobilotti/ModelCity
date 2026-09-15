# Getting started

This guide takes an ordinary pandas table through normalization, analysis, and
plotting. It does not depend on either example dataset included in the
repository.

## Requirements and installation

ModelCity requires Python 3.9 or newer. Create an isolated environment from the
repository root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

The optional `geo` extra installs `pyproj`, which is needed only by readers that
convert projected coordinates:

```bash
python -m pip install -e ".[dev,geo]"
```

Check the installation:

```bash
python -c "import modelcity; print(modelcity.__version__)"
pytest
```

## Build a settlement table

ModelCity accepts a pandas `DataFrame`. At minimum, each row needs:

- a settlement identifier;
- a non-negative area or population value;
- a start or observation year.

An end year distinguishes interval records from snapshots.

```python
import pandas as pd

frame = pd.DataFrame(
    {
        "site_code": ["A", "A", "B", "B", "C", "C"],
        "occupied_ha": [8.0, 15.0, 5.0, 9.0, 0.0, 22.0],
        "period_start_bp": [3200, 2900, 3200, 2900, 3200, 2900],
        "period_end_bp": [2900, 2600, 2900, 2600, 2900, 2600],
        "region": ["north", "north", "south", "south", "south", "south"],
    }
)
```

A size of zero means that the settlement was not occupied in that period.
Negative source sizes are rejected. Always use a positive `size_scale`.
The example uses shared period bounds, as required by period-level summaries.

## Normalize the table

Map source columns to semantic roles rather than renaming the input:

```python
import modelcity as mc

settlements = mc.as_settlements(
    frame,
    id="site_code",
    size="occupied_ha",
    start="period_start_bp",
    end="period_end_bp",
    group="region",
    size_type="area",
    time_scale="BP",
)

print(settlements)
print(settlements.data.head())
```

The resulting `Settlements` object stores a canonical table plus the size,
time, and original role metadata needed by every downstream function.

## Calculate metrics

Metrics return pandas DataFrames:

```python
durations = mc.settlement_duration(settlements)
totals = mc.size_totals(settlements)
summary = mc.size_summary(settlements)

print(durations[["id", "duration", "max_size"]])
print(totals[["label", "total_size", "n_settlements"]])
```

This separation makes it possible to inspect or export the analytical values
without reconstructing them from a chart.

## Create and save a figure

Plot functions return ordinary Matplotlib figures:

```python
figure = mc.plot_size_summary(settlements)
mc.save_figure(figure, "outputs/size-summary.png", dpi=200)
```

Close figures explicitly in batch jobs:

```python
import matplotlib.pyplot as plt

plt.close(figure)
```

## Change date labels

ModelCity stores years internally in astronomical numbering. Changing display
labels cannot alter numerical results or stored years, although display `label`
columns in newly returned tables follow the active setting:

```python
with mc.year_display("BCAD"):
    figure = mc.plot_events(settlements)
```

Use `mc.set_year_display("BCAD")` to change labels globally and
`mc.set_year_display("BP")` to restore the default.

## Next steps

- Understand the [canonical data model](data-model.md).
- Follow the complete [analysis workflow](workflow.md).
- Learn how to [adapt another dataset](adapting-datasets.md).
- Choose from the available [metrics and plots](metrics-and-plots.md).
