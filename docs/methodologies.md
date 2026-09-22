# Methodologies by stage

ModelCity splits work into stages. Source-specific loading stays outside the
package; everything below is the shared library path after you have a pandas
DataFrame.

```text
Data preparation  →  Analysis  →  Visualisation  →  Output
     (normalize)      (tables)      (figures)       (save / labels)
```

For formulas, caveats, and bandwidth guidance, see
[Metrics and plots](metrics-and-plots.md) and the
[Canonical data model](data-model.md). This page is the inventory of public
functions and the stage each belongs to.

## Stage overview

| Stage | Role | Typical inputs | Typical outputs |
| --- | --- | --- | --- |
| Data preparation | Map columns, scale size, normalise time, build a canonical table | `DataFrame` + role mapping | `Settlements` |
| Time helpers | Convert and label years; infer snapshot interval bounds | year values / display choice | arrays, strings, bounds |
| Analysis | Compute summaries used for comparison across datasets | `Settlements` | pandas `DataFrame` |
| Visualisation | Render analysis results as figures | `Settlements` | Matplotlib `Figure` |
| Output | Persist figures and control axis labels | figure / display mode | file path / label text |

---

## Data preparation

These functions turn a source table into one canonical `Settlements` object.
They do not compute research summaries.

| Function / type | What it does |
| --- | --- |
| `as_settlements(frame, ...)` | Resolves column roles, coerces numerics, applies size scale and time scale, builds interval or snapshot records, and returns a `Settlements` object. |
| `Settlements` | Canonical container: normalised table plus `SizeSpec` and `TimeSpec` metadata used by metrics and plots. |
| `SizeSpec` | Declares whether size is area or population, the unit, the scale factor already applied, and default urban thresholds. |
| `TimeSpec` | Declares input time scale (`BP` or `CE`), record type (`interval` or `snapshot`), and whether years were converted to astronomical form. |
| `ColumnRoles` | Holds the resolved mapping from canonical roles (`id`, `size`, `start`, …) to source column names. |
| `SchemaError` | Raised when required roles cannot be resolved or values fail validation. |
| `REQUIRED_ROLES` | Canonical roles that must be mapped: `id`, `size`, `start`. |
| `OPTIONAL_ROLES` | Canonical roles that may be omitted: `end`, `name`, `group`. |
| `ALL_ROLES` | Union of required and optional roles. |
| `SIZE_TYPES` | Allowed size semantics: `area` and `population`. |

**Stage note.** Calling `as_settlements()` is the only required library step before
analysis. Dataset adapters (for example in `scripts/run_unified.py`) belong
upstream of this stage.

---

## Time helpers

Supporting utilities for dates. Some run inside data preparation (`to_astronomical`,
snapshot bounds); others only affect display labels.

| Function | What it does |
| --- | --- |
| `to_astronomical(bp)` | Converts years BP to astronomical calendar years (`year = 1950 − bp`) for internal arithmetic. |
| `to_bp(year)` | Inverse of `to_astronomical`. |
| `bcad_to_bp(x, ...)` | Historical BC/AD ↔ BP conversion matching the project’s R helper (no year zero; R-compatible post-1950 quirk). |
| `ce_to_bp(ce)` | Alias path through `bcad_to_bp` for calendar-year inputs. |
| `bp_to_ce(bp)` | Inverse historical conversion from BP to BC/AD labels. |
| `format_year(value, ...)` | Formats one year for display as BP or BC/AD without changing stored values. |
| `format_period(start, end, ...)` | Formats a period label, collapsing a shared era suffix when possible. |
| `get_year_display()` / `set_year_display(display)` | Reads or sets the global label mode (`BP` default, or `BCAD`). |
| `year_display(display)` | Context manager that temporarily switches label mode for one block of plots. |
| `YEAR_DISPLAYS` | Allowed display modes. |

Snapshot interval inference (`snapshot_bounds`) runs inside `as_settlements`
when no end column is supplied: observation years become contiguous bins via
midpoints between sorted years, with outer half-gap extension.

---

## Analysis

Metric functions return pandas DataFrames. All settlement- and period-level
metrics keep only rows with `size > 0` (zero means unoccupied).

| Function | What it does |
| --- | --- |
| `settlement_duration(settlements)` | Per settlement: first occupied start, last occupied end, duration, max size, and related fields. |
| `settlement_events(settlements, bw=150, n=512)` | Kernel densities of foundation years and abandonment years on a shared grid. |
| `settlement_churn(settlements, bw=150, n=512)` | From event densities: `churn = foundations + abandonments` and `net = foundations − abandonments`. |
| `size_totals(settlements)` | Per period: total occupied size and count of unique occupied settlements. |
| `size_summary(settlements)` | Per period: median, mean, min, max, count, and 10th/90th percentiles of occupied sizes. |
| `rank_size(settlements)` | Per period: settlements ranked by size descending, plus Zipf ideal `largest / rank`. |
| `urban_share(settlements, thresholds=None)` | Per period: absolute and percent share of size at or above large-centre thresholds. |
| `kde(values, bw, n=512, ...)` | Low-level Gaussian KDE with absolute bandwidth (R `density`-style); used by event and duration analyses. |

**Stage note.** Analysis does not draw figures. Pair each table with the plot
of the same name when you need graphics.

Default urban thresholds come from `SizeSpec`: 20 and 40 ha for area data, or
1000 and 2000 inhabitants for population data.

---

## Visualisation

Plot functions call the matching analysis routine and return a Matplotlib
figure. They do not change the underlying numbers except for presentation
choices (titles, dual axes, optional log scale, peak markers).

| Function | What it does | Analysis it uses |
| --- | --- | --- |
| `plot_duration(...)` | Density of settlement durations; optional peak marks; median in the title. | `settlement_duration` + `kde` |
| `plot_events(...)` | Foundation and abandonment density curves over time. | `settlement_events` |
| `plot_churn(...)` | Churn and net density curves over time. | `settlement_churn` |
| `plot_size_totals(...)` | Total system size by period and settlement count on a second axis. | `size_totals` |
| `plot_size_summary(...)` | Median, mean, percentile band, and maximum size by period. | `size_summary` |
| `plot_rank_size(...)` | Rank–size panels per period with Zipf reference line. | `rank_size` |
| `plot_urban_share(...)` | Large-centre absolute size and percent share by period. | `urban_share` |

**Stage note.** Wording (“site” vs “city”) follows `size_type`. Year axis labels
follow the current year-display setting; duration plots use elapsed years and
are unaffected by BP/BC-AD mode.

---

## Output

| Function | What it does |
| --- | --- |
| `save_figure(fig, path, dpi=150)` | Writes a figure to disk (creating parent directories as needed) and returns the path. |

Use `set_year_display` or `year_display` before or around plotting when you want
BC/AD labels on time axes. That is presentation only; stored years stay
astronomical where conversion applied.

---

## Recommended order of use

1. **Prepare:** map your table with `as_settlements(...)`.
2. **Inspect:** print the `Settlements` object; check `size_type`, periods, and record type.
3. **Analyse:** call the metric functions you need; keep the DataFrames for tables or checks.
4. **Visualise:** call the matching `plot_*` functions.
5. **Output:** `save_figure(...)`; optionally export metric tables with `DataFrame.to_csv`.

The demonstration runner [`scripts/run_unified.py`](../scripts/run_unified.py)
follows that order after its script-local adapters produce `Settlements`
objects.
