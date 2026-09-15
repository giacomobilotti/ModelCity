# Adapting a new dataset

ModelCity does not require a fixed CSV schema. Adaptation means describing
which source columns play the canonical roles and declaring the meaning of size
and time.

## Preparation checklist

Before calling `as_settlements()`, determine:

1. What uniquely identifies a settlement?
2. Does each row describe an interval or a dated snapshot?
3. Is size an area or a population?
4. What unit and scale factor does the size column use?
5. Are dates BP or calendar-year values?
6. Does zero mean unoccupied, and how are unknown values encoded?
7. Can a settlement appear more than once in the same period?

Inspect the source rather than relying on column names:

```python
import pandas as pd

frame = pd.read_csv("data/new-dataset.csv")
print(frame.dtypes)
print(frame.head())
print(frame.isna().sum())
```

## Recipe 1: renamed interval columns

```python
import modelcity as mc

settlements = mc.as_settlements(
    frame,
    id="catalogue_number",
    name="place_label",
    size="surface_hectares",
    start="earliest_bp",
    end="latest_bp",
    group="survey_zone",
    size_type="area",
    time_scale="BP",
)
```

The source columns are preserved through `settlements.roles.mapping()`, while
the canonical frame uses `id`, `name`, `group`, `t_start`, `t_end`, and `size`.

## Recipe 2: population stored in thousands

```python
settlements = mc.as_settlements(
    frame,
    id="city_code",
    size="population_000s",
    start="census_year",
    size_type="population",
    size_scale=1000,
    size_unit="inhabitants",
    time_scale="CE",
)
```

`size_scale` is applied once during construction. Thresholds and outputs then
use inhabitants, not thousands. Always use a positive multiplier; the current
constructor does not reject a negative `size_scale`.

## Recipe 3: missing optional columns

Only `id`, `size`, and `start` are required. Omit `name`, `group`, and `end`
when unavailable:

```python
settlements = mc.as_settlements(
    frame,
    id="site",
    size="area",
    start="year",
    size_type="area",
    time_scale="CE",
)
```

ModelCity copies `id` into `name`, fills `group` with missing values, and
interprets the rows as snapshots because no `end` was mapped.

## Recipe 4: canonical column names

When the source already contains `id`, `size`, `start`, and optionally `end`,
only semantic metadata is needed:

```python
settlements = mc.as_settlements(
    frame,
    size_type="population",
    time_scale="CE",
)
```

Canonical role names are resolved automatically, including case-only
differences. Other aliases are suggested in errors but never silently selected.

## Recipe 5: preserve extra metadata

Columns that are not used as roles are copied into the canonical frame unless
their names collide with canonical columns:

```python
settlements = mc.as_settlements(
    frame,
    id="site_id",
    size="area_ha",
    start="start_bp",
    end="end_bp",
    size_type="area",
    time_scale="BP",
)

coordinates = settlements.data[["id", "longitude", "latitude"]]
```

This is useful for later joins, faceting, mapping, or provenance. Core metrics
ignore these columns unless a function explicitly documents otherwise.
Rename a source column named `t_start`, `t_end`, `id`, `name`, `group`, or
`size` before construction if it is not serving that canonical role and must be
preserved.

## Recipe 6: unsupported time conventions

ModelCity directly accepts BP and calendar-year values. Convert other
conventions before normalization:

```python
frame["year_ce"] = convert_source_dates(frame["source_date"])

settlements = mc.as_settlements(
    frame,
    id="site",
    size="population",
    start="year_ce",
    size_type="population",
    time_scale="CE",
)
```

Document the conversion and retain the original date column. For historical BC
input, remember that ModelCity's internal arithmetic uses astronomical years
with a year zero. The helpers `bcad_to_bp()` and `to_astronomical()` can be
combined when an explicit conversion is needed. The compatibility
`bcad_to_bp()`/`ce_to_bp()` implementation reproduces an original R behavior
that gives incorrect positive BP values for dates after 1950 CE; convert modern
dates directly instead.

## Cleaning missing and duplicate values

`as_settlements()` can coerce numeric size and time columns and, by default,
drops rows missing required numeric values. Clean missing identifiers before
construction because identifiers are represented as strings in the canonical
table:

```python
frame = frame.dropna(subset=["site_id"])
```

Decide how duplicates should be resolved. Repeated rows affect size sums,
distributions, duration row counts, rank-size entries, and threshold shares;
only settlement-count fields use unique identifiers. Use one row per settlement
per exact period unless duplication is intentional. Aggregate accidental
duplicates before normalization:

```python
frame = (
    frame.groupby(["site_id", "year"], as_index=False)
    .agg(population=("population", "max"))
)
```

Use a different aggregation only when justified by the source.

## Diagnosing a mapping error

An invalid mapping raises `SchemaError` with all unresolved roles:

```python
try:
    settlements = mc.as_settlements(
        frame,
        id="missing_id",
        size="pop",
        start="date",
        size_type="population",
    )
except mc.SchemaError as error:
    print(error)
```

The message reports:

- what each unresolved role means;
- which column was requested;
- all available columns;
- a likely alternative where one can be identified;
- a copyable `as_settlements(...)` repair template.

## Write a reusable reader

When a schema will be used repeatedly, isolate its mapping in a thin reader:

```python
def read_my_dataset(path, **overrides):
    frame = pd.read_csv(path)
    options = {
        "id": "site_code",
        "size": "area_ha",
        "start": "start_bp",
        "end": "end_bp",
        "size_type": "area",
        "time_scale": "BP",
    }
    options.update(overrides)
    return mc.as_settlements(frame, **options)
```

The reader should perform source-specific repair and metadata declaration only.
Keep metrics and plotting in shared code.

## Included readers as examples

[`read_yautepec()`](../src/modelcity/readers.py) maps an interval table with
areas in hectares and BP dates.

[`read_viabundus()`](../src/modelcity/readers.py) demonstrates a more involved
snapshot reader: it repairs an R-exported CSV, parses coordinates, optionally
converts them with `pyproj`, scales populations from thousands, and then calls
the same `as_settlements()` constructor.

These readers illustrate adapter patterns; neither schema is privileged by the
core package.

After adapting a source, follow the [end-to-end workflow](workflow.md) and the
[reproducibility checklist](reproducibility.md).
