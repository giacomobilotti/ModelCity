# Python API

This page summarizes the supported top-level API. Import these names from
`modelcity`; submodule helpers not re-exported by
[`src/modelcity/__init__.py`](../src/modelcity/__init__.py) should be treated as
internal unless documented otherwise.

Docstrings in the source remain the parameter-level source of truth:

```python
import modelcity as mc

help(mc.as_settlements)
help(mc.settlement_churn)
```

## Construction and validation

### `as_settlements(data, *, ...)`

Normalize a pandas DataFrame into a `Settlements` object.

Important keyword arguments:

- role mappings: `id`, `size`, `start`, `end`, `name`, `group`;
- `size_type`: `"area"` or `"population"` (required);
- `time_scale`: `"BP"` or `"CE"` (default `"CE"`);
- `size_scale`: numeric multiplier (default `1.0`);
- `size_unit`: unit-label override;
- `drop_missing`: drop rows missing required numeric normalized values (default
  `True`). Clean missing identifiers before construction because IDs are
  converted to strings.

See [Canonical data model](data-model.md) for semantics.

### `Settlements`

Wrapper around the canonical DataFrame and its analytical metadata.

Core attributes:

- `data`: canonical pandas DataFrame;
- `size`: `SizeSpec`;
- `time`: `TimeSpec`;
- `roles`: `ColumnRoles`.

Properties and methods:

- `n_settlements`;
- `n_periods`;
- `is_snapshot`;
- `size_label`;
- `has_group`;
- `periods()`;
- `occupied()`;
- `years_bp()`;
- `with_data(data)`.

### Metadata classes

- `SizeSpec(size_type, unit, scale=1.0)`: exposes `is_area`, `label`, and
  `default_thresholds`.
- `TimeSpec(time_scale, record_type)`: exposes `is_snapshot`.
- `ColumnRoles(id, size, start, end=None, name=None, group=None)`: exposes
  `mapping()` and `has_end`.
- `SchemaError`: raised for invalid mappings or unusable normalized values.

### Constants

- `SIZE_TYPES`: supported size semantics;
- `REQUIRED_ROLES`: `id`, `size`, `start`;
- `OPTIONAL_ROLES`: `end`, `name`, `group`;
- `ALL_ROLES`: required and optional roles together.

Source-specific adapters for the repository example CSVs live in
[`scripts/run_unified.py`](../scripts/run_unified.py), not in the package API.
New datasets should use `as_settlements()` directly or define a similarly thin
adapter outside the library.

## Metrics

All metric functions return pandas DataFrames.

- `settlement_duration(settlements)`;
- `settlement_events(settlements, bw=150.0, n=512)`;
- `settlement_churn(settlements, bw=150.0, n=512)`;
- `size_totals(settlements)`;
- `size_summary(settlements)`;
- `rank_size(settlements)`;
- `urban_share(settlements, thresholds=None)`;
- `kde(values, bw, n=512, cut=3.0, lower=None, upper=None)`.

See [Metrics and plots](metrics-and-plots.md) for outputs and interpretation.
With no occupied rows, `urban_share()` currently returns a columnless frame and
`plot_urban_share()` cannot render it.

## Plots

All plot functions accept a `Settlements` object as their first argument and
return a `matplotlib.figure.Figure`.

- `plot_duration(settlements, bw=150.0, mark_peaks=True, figsize=(10, 7.5))`;
- `plot_events(settlements, bw=150.0, figsize=(10, 7.5))`;
- `plot_churn(settlements, bw=150.0, figsize=(10, 7.5))`;
- `plot_size_totals(settlements, figsize=(10, 7.5))`;
- `plot_size_summary(settlements, figsize=(10, 7.5), log_scale=False)`;
- `plot_rank_size(settlements, ncol=4, figsize=None)`;
- `plot_urban_share(settlements, thresholds=None, figsize=(10, 7.5))`.

### `save_figure(fig, path, dpi=150)`

Save a Matplotlib figure, creating parent directories as needed. Returns the
destination as a string.

## Time conversion

### Arithmetic scales

- `to_astronomical(bp)`: BP to astronomical calendar years;
- `to_bp(year)`: astronomical calendar years to BP.

Use these linear conversions for calculations.

### Historical BC/AD compatibility

- `bcad_to_bp(x, bc_to_bp=True)`: conversion compatible with the project's R
  implementation;
- `ce_to_bp(ce)`: historical BC/AD to BP;
- `bp_to_ce(bp)`: BP to historical BC/AD values for labels.

Historical BC/AD has no year zero. Astronomical numbering does.
These helpers intentionally match the project's original R conversion,
including its incorrect positive-BP result for post-1950 CE dates. Do not use
`bcad_to_bp()` or `ce_to_bp()` for those modern dates.

### Formatting

- `format_year(year, suffix=True, display=None)`;
- `format_period(start, end, display=None)`.

Both accept astronomical years and use the active display unless explicitly
overridden.

## Display configuration

- `YEAR_DISPLAYS`: `("BP", "BCAD")`;
- `get_year_display()`: current global display;
- `set_year_display(display)`: set display and return the previous value;
- `year_display(display)`: context manager that restores the previous value.

```python
with mc.year_display("BCAD"):
    figure = mc.plot_events(settlements)
```

Display configuration affects labels only.

## Version

The package version is available as:

```python
mc.__version__
```

The current package version is `0.1.0`.
