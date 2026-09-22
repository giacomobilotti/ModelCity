# End-to-end workflow

The reusable ModelCity workflow has one source-specific step—loading and
mapping a table—and a shared analytical path after normalization.

```text
load source data
    ├── interval records: map start + end
    └── snapshot records: map observation year only
             ↓
       as_settlements()          ← data preparation
             ↓
 inspect canonical data and assumptions
             ↓
 calculate DataFrame metrics     ← analysis
             ↓
 create and save Matplotlib figures  ← visualisation / output
```

A function-by-function map of these stages is in
[Methodologies by stage](methodologies.md).

## 1. Load a DataFrame

Use any reader that produces a pandas DataFrame:

```python
from pathlib import Path
import pandas as pd

source = Path("data/my-settlements.csv")
frame = pd.read_csv(source)
print(frame.columns.tolist())
print(frame.head())
```

Clean source-specific encodings before normalization. In particular, decide how
missing sizes, duplicate observations, uncertainty, and zero values should be
represented. ModelCity validates its required roles but does not infer the
research meaning of unusual source values.

## 2. Declare mappings and metadata

### Interval route

Use intervals when rows state the period over which a value applies:

```python
import modelcity as mc

settlements = mc.as_settlements(
    frame,
    id="settlement_id",
    name="settlement_name",
    size="occupied_area",
    start="period_start",
    end="period_end",
    group="region",
    size_type="area",
    size_unit="km²",
    time_scale="CE",
)
```

Period metrics group exact start/end pairs. Harmonize intervals into shared
period bins before normalization; overlapping records with different bounds
are not automatically aligned.

### Snapshot route

Use snapshots when rows state a value observed at one date:

```python
settlements = mc.as_settlements(
    frame,
    id="settlement_id",
    name="settlement_name",
    size="population_estimate",
    start="observation_year",
    group="region",
    size_type="population",
    time_scale="CE",
)
```

With no `end` mapping, ModelCity infers period bounds halfway between distinct
observation years. Review that assumption before interpreting durations,
foundations, or abandonments.

## 3. Inspect the normalized object

Do not proceed directly from construction to publication. Check the normalized
shape and metadata:

```python
print(settlements)
print(settlements.roles.mapping())
print(settlements.size)
print(settlements.time)
print(settlements.periods())
print(settlements.data.head())
```

Useful checks:

```python
assert settlements.data["size"].ge(0).all()
assert settlements.data["t_start"].le(settlements.data["t_end"]).all()
assert settlements.n_settlements > 0
```

Confirm that:

- `size_scale` has produced the intended real-world unit;
- `size_scale` is positive;
- source zeros truly mean absence rather than missing data;
- interval bounds are chronological after BP conversion;
- inferred snapshot periods are appropriate;
- the number of unique settlements and periods is plausible.

`group` is descriptive metadata in the current API. Metrics do not facet by it,
and duration output retains only the first group encountered per settlement.

## 4. Calculate metrics first

All analytical functions return DataFrames:

```python
durations = mc.settlement_duration(settlements)
events = mc.settlement_events(settlements, bw=100)
churn = mc.settlement_churn(settlements, bw=100)
totals = mc.size_totals(settlements)
summary = mc.size_summary(settlements)
ranked = mc.rank_size(settlements)
shares = mc.urban_share(settlements, thresholds=[10, 25])
```

Inspect and export these values separately from plots:

```python
output = Path("outputs")
output.mkdir(parents=True, exist_ok=True)

durations.to_csv(output / "durations.csv", index=False)
totals.to_csv(output / "size-totals.csv", index=False)
shares.to_csv(output / "large-centre-share.csv", index=False)
```

Choose kernel bandwidths and large-centre thresholds based on the temporal
resolution, units, and research question. Defaults provide comparable behavior;
they are not universal substantive definitions.

## 5. Create figures

Each plot function accepts the same `Settlements` object and returns a
Matplotlib `Figure`:

```python
plots = {
    "duration": mc.plot_duration(settlements, bw=100),
    "events": mc.plot_events(settlements, bw=100),
    "churn": mc.plot_churn(settlements, bw=100),
    "size-totals": mc.plot_size_totals(settlements),
    "size-summary": mc.plot_size_summary(settlements),
    "rank-size": mc.plot_rank_size(settlements),
    "large-centre-share": mc.plot_urban_share(
        settlements, thresholds=[10, 25]
    ),
}
```

Save and close figures in batch work:

```python
import matplotlib.pyplot as plt

for name, figure in plots.items():
    mc.save_figure(figure, output / f"{name}.png", dpi=200)
    plt.close(figure)
```

## 6. Choose display labels

BP is the default display. Switch labels for one operation:

```python
with mc.year_display("BCAD"):
    figure = mc.plot_size_totals(settlements)
```

The display setting changes tick and period labels only. Dates remain
astronomical internally, and calculated values do not move.

## 7. Turn the workflow into a project script

Keep loading separate from shared analysis:

```python
def load_settlements(path):
    frame = pd.read_csv(path)
    return mc.as_settlements(
        frame,
        id="settlement_id",
        size="population",
        start="year",
        size_type="population",
        time_scale="CE",
    )


def analyze(settlements, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    mc.size_totals(settlements).to_csv(
        output / "size-totals.csv", index=False
    )
    figure = mc.plot_size_totals(settlements)
    mc.save_figure(figure, output / "size-totals.png")
```

This boundary makes adding another dataset a small loader change rather than a
copy of the analytical code. The included
[`scripts/run_unified.py`](../scripts/run_unified.py) keeps those loaders
outside the package and demonstrates the pattern with two adapters and one
shared set of metrics and plots.

## Assumptions to report

For a reproducible analysis, record:

- source data version and preprocessing;
- role mappings and omitted records;
- interval versus snapshot interpretation;
- source and display time scales;
- size type, unit, and scale factor;
- bandwidths and urban thresholds;
- ModelCity version and Python environment.

See [Reproducibility](reproducibility.md) for a checklist.
