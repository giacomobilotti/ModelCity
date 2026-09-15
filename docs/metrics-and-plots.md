# Metrics and plots

ModelCity keeps numerical summaries separate from visualization. Metric
functions return pandas DataFrames; their paired plot functions return
Matplotlib figures.

All settlement-level and period-level metrics consider only rows where
`size > 0`. A zero therefore means unoccupied, not an observed settlement of
zero extent or population.

## Persistence

### `settlement_duration(settlements)`

Question: how long did each settlement persist?

For each identifier, the function spans the earliest occupied `t_start` to the
latest occupied `t_end`. It returns:

- `id`, `name`, and `group`;
- `first_year` and `last_year`;
- `max_size`;
- `n_periods`, which currently counts occupied rows rather than distinct period
  bounds;
- `duration`.

For snapshots, the first and last years are the inferred outer period bounds,
not merely the observation dates.

### `plot_duration(settlements, bw=150, mark_peaks=True, figsize=...)`

Plots a kernel density of durations and reports the median in the title.
`mark_peaks` annotates local density maxima. The x axis is elapsed years and is
not affected by the BP/BC-AD display setting.

## Foundations, abandonments, and turnover

### `settlement_events(settlements, bw=150, n=512)`

Question: when did settlements first appear and last disappear?

The result contains `year`, `density`, and `event`, with events labeled
`foundations` or `abandonments`. Both curves use the same grid and bandwidth.

### `settlement_churn(settlements, bw=150, n=512)`

Question: how did the relative timing density of first and last appearances
compare?

The result contains:

- `foundations` and `abandonments`: event densities;
- `churn`: the sum of the two independently normalized densities;
- `net`: foundation density minus abandonment density.

### `plot_events(...)` and `plot_churn(...)`

`plot_events` compares the two event-density curves. `plot_churn` overlays their
sum and difference. These values are density comparisons, not event counts.
They use only each settlement's first and last occupied bounds, so intermediate
disappearance and reappearance are not represented. Positive `net` means first
appearances have greater relative density at that date; it does not by itself
measure population growth or a change in the number of active settlements.

Bandwidth is absolute in years. A larger bandwidth emphasizes broad trends; a
smaller value exposes shorter fluctuations but may overstate noise. Choose it
in relation to dating precision and observation spacing, and report it.

## System size and settlement counts

### `size_totals(settlements)`

Question: how large was the complete occupied settlement system in each period?

It returns period bounds and labels plus:

- `total_size`: sum of occupied areas or populations;
- `n_settlements`: number of unique occupied identifiers.

Periods with no occupied rows are retained with zeros.

Period summaries group exact `(t_start, t_end)` pairs. They require harmonized
shared period bins for a system-wide interpretation; arbitrary overlapping
intervals are not aligned automatically.

### `plot_size_totals(settlements, figsize=...)`

Shows total size as period segments and settlement count on a second y axis.
Units and wording come from `SizeSpec`.

## Size distribution

### `size_summary(settlements)`

Question: what did the distribution of occupied settlement sizes look like in
each period?

It returns `median`, `mean`, `min`, `max`, `n`, `q10`, and `q90`, together with
period bounds, midpoint, and label.

### `plot_size_summary(settlements, figsize=..., log_scale=False)`

Shows median and mean, the 10th–90th percentile band, and the largest
settlement. Set `log_scale=True` for highly skewed positive sizes. Zero-size
rows are already excluded.

## Rank-size structure

### `rank_size(settlements)`

Question: how hierarchical was the settlement-size distribution?

Within each period, occupied settlements are sorted from largest to smallest.
The result adds:

- `rank`;
- `ideal`, the Zipf expectation `largest_size / rank`;
- period midpoint and display label.

Repeated identifiers in the same period are not aggregated before ranking.
Inputs should therefore have one row per settlement-period unless repeated rows
are substantively intended. Duplicates also affect size totals, distributions,
duration row counts, and threshold shares; only settlement-count fields use
unique identifiers.

### `plot_rank_size(settlements, ncol=4, figsize=None)`

Creates one log-log panel per period and compares observed sizes with the Zipf
expectation. It raises `ValueError` when no periods are available. `ncol`
controls the panel layout.

## Large-centre share

### `urban_share(settlements, thresholds=None)`

Question: how much of the system's total size was concentrated in settlements
at or above selected thresholds?

For every threshold `x`, the result includes:

- `urban_x`: combined size at or above the threshold;
- `n_urban_x`: number of qualifying rows;
- `percent_x`: percentage of total occupied size.

It also includes total size, settlement count, and period columns. Custom
thresholds are strongly recommended when the package defaults do not match the
research definition:

```python
shares = mc.urban_share(settlements, thresholds=[5000, 10000, 25000])
```

Defaults are 20 and 40 hectares for area data, and 1,000 and 2,000 inhabitants
for population data. These values provide package behavior; they are not a
universal definition of urbanism.

At least one occupied row is required for `plot_urban_share()`. An all-zero or
empty dataset currently produces a columnless metric frame that cannot be
plotted.

### `plot_urban_share(settlements, thresholds=None, figsize=...)`

Plots absolute size above each threshold on the primary axis and percentage of
the system total as dashed lines on the secondary axis.

Despite the function name, its general interpretation is **large-centre
share**. Whether a threshold represents “urban” status is a domain decision
that should be justified by the analysis.

Current plot wording also calls area-based settlements “sites” and
population-based settlements “cities.” This is a fixed presentation convention
in the plotting module.

## Standalone density estimation

### `kde(values, bw, n=512, cut=3, lower=None, upper=None)`

Returns `year` and `density` columns for a Gaussian kernel-density estimate.
`bw` is an absolute kernel bandwidth, matching the interpretation of R's
`density(..., bw=...)` rather than SciPy's default relative bandwidth.

Unless explicit bounds are supplied, the grid extends `cut × bw` beyond the
data. Empty input returns an empty frame; a single value or constant series is
handled with a directly evaluated Gaussian kernel.

## Saving and customizing figures

Every plot function returns a regular `matplotlib.figure.Figure`:

```python
figure = mc.plot_rank_size(settlements, ncol=3)
mc.save_figure(figure, "outputs/rank-size.png", dpi=200)
```

You may edit titles, labels, or artists through standard Matplotlib APIs before
saving. `save_figure()` creates parent directories and writes with a white
background and tight bounding box.

See the [Python API](python-api.md) for the stable exported function list.
