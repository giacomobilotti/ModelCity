# Reproducibility

A reproducible ModelCity analysis records the source transformation and
analytical assumptions, not only the plotting command.

## Record the environment

Create an isolated environment and capture exact installed versions:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
python --version
python -m pip freeze > environment.txt
```

For a published analysis, also record:

- ModelCity version or Git commit;
- operating system;
- input-data version, checksum, or persistent identifier;
- the command used to run the analysis.

The project metadata defines lower dependency bounds, not a fully locked
environment. An exported environment or external lock file is therefore needed
for exact recreation.

## Preserve source and normalized assumptions

Store the `as_settlements()` configuration next to the analysis:

```python
MAPPING = {
    "id": "site_code",
    "size": "area_ha",
    "start": "start_bp",
    "end": "end_bp",
    "group": "region",
    "size_type": "area",
    "time_scale": "BP",
}

settlements = mc.as_settlements(frame, **MAPPING)
```

Report:

- interval or snapshot interpretation;
- handling of missing, duplicate, uncertain, and zero values;
- size type, source unit, `size_scale`, and output unit;
- source date convention and any preprocessing conversion;
- grouping semantics;
- excluded rows or periods.

For snapshots, state that ModelCity inferred interval bounds from adjacent
observation-year midpoints.

## Record analytical parameters

Defaults are still assumptions. Save the values used for:

- KDE bandwidth (`bw`) and grid size (`n`);
- large-centre thresholds;
- rank-size period definitions;
- year display;
- output DPI and file format.

Prefer named constants in the analysis script:

```python
EVENT_BANDWIDTH = 100
LARGE_CENTRE_THRESHOLDS = [20, 40]

events = mc.settlement_events(settlements, bw=EVENT_BANDWIDTH)
shares = mc.urban_share(
    settlements, thresholds=LARGE_CENTRE_THRESHOLDS
)
```

## Test before generating outputs

Run the repository suite:

```bash
pytest
```

The tests cover:

- role resolution and actionable schema errors;
- interval and snapshot normalization;
- BP, astronomical, and BC/AD conversion;
- duration, turnover, size, rank-size, and threshold metrics;
- all seven plot functions for both area and population data;
- integration with the included example datasets when they are available;
- parity of KDE behavior and key time conversions with the original R work.

Passing tests demonstrates that the package behaves as asserted. It does not
validate the substantive suitability of a new dataset's mappings, periods, or
thresholds.

## Validate a new dataset

Add project-specific assertions before plotting:

```python
analysis_frame = frame.dropna(
    subset=["site_code", "area_ha", "start_bp", "end_bp"]
)
assert settlements.n_settlements == analysis_frame["site_code"].nunique()
assert settlements.data["size"].ge(0).all()
assert settlements.data["t_start"].le(settlements.data["t_end"]).all()
assert settlements.periods()["midpoint"].is_monotonic_increasing
```

Also compare:

- input and canonical row counts;
- the first and last source dates after conversion;
- a few manually calculated durations;
- period totals against independent sums;
- threshold classifications near boundary values.

If the dataset will become a maintained adapter, add compact fixtures and
golden-value tests under [`tests/`](../tests/).

## Separate values from figures

Export metric tables alongside images:

```python
from pathlib import Path
import matplotlib.pyplot as plt

output = Path("outputs")
output.mkdir(parents=True, exist_ok=True)

totals = mc.size_totals(settlements)
totals.to_csv(output / "size-totals.csv", index=False)

figure = mc.plot_size_totals(settlements)
mc.save_figure(figure, output / "size-totals.png", dpi=200)
plt.close(figure)
```

Tables make numerical comparison possible even when Matplotlib, fonts, or image
metadata differ between environments.

## Included reproducible demonstration

[`scripts/run_unified.py`](../scripts/run_unified.py) loads each included
example through a source-specific adapter (defined in that script, not in the
package) and applies one shared set of seven plots and four exported metric
tables:

```bash
python3 scripts/run_unified.py
python3 scripts/run_unified.py --dataset viabundus
```

Outputs are written under:

- `figures/yautepec/unified/`;
- `figures/viabundus/unified/`.

This script is evidence that the same API can operate on contrasting schemas.
It is not a generic command-line interface for arbitrary datasets; general
projects should create a small adapter and runner as shown in the
[workflow guide](workflow.md).

## Suggested output manifest

For each run, retain a small text or JSON manifest containing:

```text
modelcity_version
git_commit
python_version
input_path_and_checksum
mapping
size_type_scale_and_unit
time_scale_and_record_type
bandwidths_and_thresholds
year_display
run_timestamp
```

ModelCity does not currently write this manifest automatically. Treat it as a
recommended project-level practice.

## Known boundaries

- Snapshot bounds are inferred, not observed.
- Kernel bandwidth and large-centre thresholds require domain justification.
- Historical BC notation and astronomical internal years differ at the
  no-year-zero boundary.
- The package does not currently model dating uncertainty.
- Source-specific R/Quarto analyses have additional dependencies and are not
  part of the reusable Python workflow; see
  [Legacy R workflows](legacy-r-workflows.md).
