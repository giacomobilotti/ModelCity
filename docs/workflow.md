# End-to-end workflow

The reusable ModelCity workflow has one source-specific step—loading and
mapping a table—and three shared analyses after normalization.

```text
load source data
    ├── interval records: map start + end
    └── snapshot records: map observation year only
             ↓
       as_settlements()          ← data preparation
             ↓
 inspect canonical data and assumptions
             ↓
     ┌───────────────────┼────────────────────┐
 system metrics      chronological         longevity
 and plots           sampling              (from the
 (any record type)   (interval records)    persistence draws)
```

A function-by-function map of these stages is in
[Methodologies by stage](methodologies.md). The three analyses are documented
in [Metrics and plots](metrics-and-plots.md),
[Chronological sampling](chronological-sampling.md), and
[Longevity](longevity.md). Each has a runner for the included example
datasets: [`scripts/run_unified.py`](../scripts/run_unified.py),
[`scripts/run_sampling.py`](../scripts/run_sampling.py), and
[`scripts/run_longevity.py`](../scripts/run_longevity.py).

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

## 4. System metrics

This analysis summarises the settlement system directly from the canonical
table. It accepts interval and snapshot records.

The metric functions return DataFrames:

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

## 5. System-metric figures

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

## 7. Sample dates within phases

Phase bounds leave the foundation, peak, and abandonment dates free inside
each occupied interval. Chronological sampling draws those dates many times
and carries every draw through size, growth, and persistence. The functions
require interval records.

```python
ITERATIONS = 1000
SEED = 42
SAMPLE_TERMINUS = 510
PERSISTENCE_TERMINUS = 431

trajectories = mc.sample_trajectories(
    settlements,
    ITERATIONS,
    terminus=SAMPLE_TERMINUS,
    seed=SEED,
)
by_site = mc.site_summary(trajectories)
by_period = mc.period_summary(trajectories)
persistence = mc.persistence_summary(
    trajectories, terminus=PERSISTENCE_TERMINUS
)
```

`terminus` is in the dataset's own time scale. Give the sampling terminus and
the persistence terminus separately: the examples use the start of the last
phase for sampling and its end for persistence.

Figures take the trajectory frame:

```python
site_id = settlements.data["id"].iloc[0]
plots = {
    "cagr-site": mc.plot_site_growth(trajectories, site_id),
    "spaghetti-cagr": mc.plot_trajectories(
        trajectories, site_id, value="growth_cagr"
    ),
    "cagr-by-period": mc.plot_period_growth(trajectories),
    "persistence-by-site": mc.plot_persistence(
        trajectories, terminus=PERSISTENCE_TERMINUS
    ),
}
```

A single-phase settlement has one growth point, so `plot_site_growth()` skips
it. Snapshot tables need a phase reading first: each census becomes the
period `as_settlements()` inferred around it, and unoccupied rows are dropped.
[`scripts/run_sampling.py`](../scripts/run_sampling.py) does that for
Viabundus. Date methods, summary rules, and the comparison with the R samplers
are in [Chronological sampling](chronological-sampling.md).

## 8. Fit longevity to the censored draws

Longevity starts from the persistence table of a sampling run. A settlement
still occupied at the persistence terminus has a right-censored duration. The
analysis fits a survival model to each iteration and imputes a total lifespan
for the censored rows.

```python
samples = mc.mark_censoring(
    mc.persistence_samples(trajectories, terminus=PERSISTENCE_TERMINUS)
)
fits = mc.fit_persistence(samples, distribution="gengamma")
coefficients = fits[["mu", "sigma", "Q"]].mean()

figure = mc.plot_longevity(samples, coefficients=coefficients)
mc.save_figure(figure, output / "longevity.png")
```

`summarize_persistence()` pools that per-draw table by settlement. Leave
`norm_pers` unchanged when imputing if the summary should keep the share from
before imputation, which is what
[`scripts/run_longevity.py`](../scripts/run_longevity.py) does. Families,
censoring, and the R correspondence are in [Longevity](longevity.md).

## 9. Turn the workflow into a project script

Keep loading separate from the shared analyses:

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
copy of the analytical code. The included runners keep those loaders in
[`scripts/run_unified.py`](../scripts/run_unified.py) and apply one analysis
each:

| Runner | Analysis | Output |
| --- | --- | --- |
| [`scripts/run_unified.py`](../scripts/run_unified.py) | System metrics and the seven shared plots | `figures/<dataset>/unified/` |
| [`scripts/run_sampling.py`](../scripts/run_sampling.py) | Chronological sampling, growth, and persistence | `figures/<dataset>/sampling-python/` |
| [`scripts/run_longevity.py`](../scripts/run_longevity.py) | Censored lifespans and imputed longevity | `figures/<dataset>/longevity-python/` |

## Assumptions to report

For a reproducible analysis, record:

- source data version and preprocessing;
- role mappings and omitted records;
- interval versus snapshot interpretation;
- source and display time scales;
- size type, unit, and scale factor;
- bandwidths and urban thresholds;
- for sampling: iterations, seed, start/end/peak methods, and both termini;
- for longevity: the survival distribution and the imputation seed;
- ModelCity version and Python environment.

See [Reproducibility](reproducibility.md) for a checklist.
