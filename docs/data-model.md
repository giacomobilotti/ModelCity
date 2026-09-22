# Canonical data model

ModelCity separates source-specific column names from the concepts needed for
analysis. `as_settlements()` maps those concepts into a canonical
`Settlements` object.

## Column roles

Required roles:

| Role | Meaning |
| --- | --- |
| `id` | Stable settlement identifier, repeated across periods |
| `size` | Numeric occupied area or population |
| `start` | Interval start or snapshot observation year |

Optional roles:

| Role | Meaning |
| --- | --- |
| `end` | Interval end; omit it for snapshot records |
| `name` | Human-readable settlement label |
| `group` | Descriptive region, phase, culture, or another category |

Pass source column names explicitly:

```python
settlements = mc.as_settlements(
    frame,
    id="location_code",
    name="location_name",
    size="estimated_population",
    start="observation_year",
    group="study_region",
    size_type="population",
    time_scale="CE",
)
```

If a source column already has a canonical role name, its mapping may be
omitted. Resolution accepts case-only differences but does not silently guess
aliases. When a role is missing, `SchemaError` reports available columns,
suggests a likely match, and supplies a repair template.

## Canonical columns

`settlements.data` always starts with:

| Column | Meaning |
| --- | --- |
| `id` | Identifier represented as text |
| `name` | Supplied label, or `id` when no name was mapped |
| `group` | Supplied grouping, or missing values |
| `t_start` | Start bound in astronomical calendar years |
| `t_end` | End bound in astronomical calendar years |
| `size` | Size after applying `size_scale` |

It also contains `t_observed`: the source observation year for snapshots and
the interval midpoint for interval records. Unmapped source columns are
normally preserved after the canonical columns, so coordinates and other
metadata remain available for later work. A source column whose name collides
with a canonical column is not copied; rename such columns before construction
if they must be retained.

The wrapper also preserves:

- `settlements.size`: `SizeSpec(size_type, unit, scale)`;
- `settlements.time`: `TimeSpec(time_scale, record_type)`;
- `settlements.roles`: the resolved source-to-role mapping.

Use `settlements.with_data(frame)` when a transformed canonical frame should
retain this metadata. Do not replace the wrapper with a bare DataFrame before
calling metrics or plots.

## Interval and snapshot records

### Intervals

Map both `start` and `end` when rows describe explicit periods:

```python
intervals = mc.as_settlements(
    frame,
    id="site",
    size="area",
    start="start_bp",
    end="end_bp",
    size_type="area",
    time_scale="BP",
)
```

For BP data, the source start is normally the larger number because BP counts
backwards. After conversion, `t_start` must still be chronologically earlier
than `t_end`; reversed records raise `SchemaError`.

Period-level metrics group records by exact `(t_start, t_end)` pairs. Inputs
should therefore use harmonized shared period bins across settlements.
Arbitrary overlapping intervals remain separate groups; ModelCity does not
split or align them automatically.

### Snapshots

Omit `end` when each row is an observation at one date:

```python
snapshots = mc.as_settlements(
    frame,
    id="city",
    size="population",
    start="year",
    size_type="population",
    time_scale="CE",
)
```

ModelCity turns distinct observation years into contiguous intervals. Interior
bounds fall halfway between adjacent observations; the first and last periods
extend outward by half their nearest gap. A dataset with only one observation
year receives a zero-width period at that year.

These inferred bounds affect duration and period-based metrics. They represent
an explicit analytical assumption, not evidence that an observation remained
constant throughout the inferred interval.

## Size metadata

`size_type` is required and must be:

- `"area"`: default unit `ha`, default large-centre thresholds 20 and 40;
- `"population"`: default unit `inhabitants`, default thresholds 1,000 and
  2,000.

Override the unit label with `size_unit`. Use `size_scale` when source values
are stored in scaled units:

```python
settlements = mc.as_settlements(
    frame,
    id="city",
    size="population_thousands",
    start="year",
    size_type="population",
    size_scale=1000,
    size_unit="inhabitants",
)
```

Scaling occurs during normalization. All metrics, thresholds, labels, and
canonical values then use the scaled unit. Use a strictly positive
`size_scale`; the current constructor does not reject a negative multiplier
after validating source sizes.

Rows with a zero size remain in the canonical table but are excluded by
`settlements.occupied()` and by metrics that analyze occupied settlements.
Negative values are invalid.

## Time representation

Set `time_scale="BP"` when source years count backwards from 1950. Set
`time_scale="CE"` for calendar-year values. Internally, all dates use
astronomical numbering, which includes year zero and therefore supports correct
elapsed-time arithmetic.

The package distinguishes three concepts:

- **BP**: years before 1950;
- **astronomical years**: internal numeric scale, where 1 BC is year 0;
- **historical BC/AD**: display notation without a year zero.

`time_scale="CE"` values are used directly as astronomical calendar numbers.
If source dates use negative historical BC notation, convert them deliberately
instead of assuming the no-year-zero correction will happen during
normalization.

Useful helpers include `to_astronomical()`, `to_bp()`, `bcad_to_bp()`,
`bp_to_ce()`, and `ce_to_bp()`. The historical conversion helpers preserve
parity with the original R function and incorrectly turn post-1950 CE dates
into positive BP values; do not use them for modern dates after 1950.
`settlements.years_bp()` produces a copy of the canonical table with readable
BP columns for export.

## Group metadata

`group` is descriptive metadata in the current API. Core period metrics do not
aggregate or facet by it. `settlement_duration()` retains the first group value
encountered for each identifier, so a group that changes through time should
not be interpreted from that output without separate handling.

## Missing values and validation

By default, rows missing size, start, or derived end values are dropped. Clean
missing identifiers before construction: identifiers are converted to strings,
so a missing source identifier can otherwise become the literal text `"nan"`.
Pass `drop_missing=False` to retain missing numeric values, but downstream
metrics may not be meaningful until they are handled.

Validation covers:

- input must be a pandas DataFrame;
- required roles must resolve to existing columns;
- size and time columns must be numeric or coercible;
- sizes cannot be negative;
- `size_type` and `time_scale` must be supported values;
- intervals cannot end before they start after time conversion.

The constructor does not currently validate that `size_scale` is positive.

See [Adapting datasets](adapting-datasets.md) for practical recipes.
