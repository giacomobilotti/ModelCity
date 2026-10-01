# Chronological sampling

A phase only bounds when something happened. A site recorded in a 400-year
phase could have been founded, reached its recorded size, or been abandoned at
any point inside it, and growth rates computed from phase midpoints hide that
uncertainty. The sampling functions draw those dates many times within the
phase bounds and carry every draw through to growth rates and persistence, so
the spread across draws shows how much a result depends on dates the record
does not fix.

The functions live in `modelcity.sampling`, are exported from `modelcity`, and
reproduce [`scripts/r/sampling.R`](../scripts/r/sampling.R) and
[`scripts/r/helpers_sampling.R`](../scripts/r/helpers_sampling.R). They need
interval records; snapshot data is rejected.

```python
import modelcity as mc

trajectories = mc.sample_trajectories(settlements, 1000, terminus=510, seed=42)

mc.site_summary(trajectories)                     # summary_site()
mc.period_summary(trajectories)                   # summary_period()
mc.persistence_summary(trajectories, terminus=431)  # persistence()

mc.save_figure(mc.plot_site_growth(trajectories, "160"), "outputs/cagr-site-160.png")
mc.save_figure(mc.plot_trajectories(trajectories, "160"), "outputs/spaghetti-cagr-site-160.png")
```

## Correspondence with R

| R | Python |
| --- | --- |
| `run_mc(sites, iterations, ...)` | `sample_trajectories(settlements, iterations, ...)` |
| `sample_dates()`, `growth_rate()` | run inside `sample_trajectories()` |
| `exp_sampling()`, `norm_peak()` | internal samplers |
| `summary_site()` | `site_summary()` |
| `summary_period()` | `period_summary()` |
| `mc_persistence(sites, terminus)` | `persistence_samples(trajectories, terminus)` |
| `persistence(sites, terminus)` | `persistence_summary(trajectories, terminus)` |
| `plot_site_cagr()` | `plot_site_growth()` |
| `plot_spaghetti()` | `plot_trajectories()` |

Arguments keep their R names and defaults: `start_method`, `end_method` and
`peak_method` (`"unif"`, `"exp"`, `"norm"`; `"uniform"`, `"exponential"` and
`"normal"` are accepted too), `perc=5`, `peak_sd_perc=0.1` and `terminus=1`
(1 BP when left as `None`). `terminus` is given in the dataset's own time
scale. Persistence has no default `terminus`; R's default of 431 BP is specific
to Yautepec.

## What is reproduced

The rules are applied in years BP and every date is floored to a whole year,
as in R:

- **Foundation**: uniform within the first phase, or the phase start minus a
  floored exponential delay below the phase length (`exp`).
- **Abandonment of a single-phase settlement**: uniform between the phase end
  and the foundation, or the foundation minus a floored exponential delay
  (`exp`), so it tends to follow the foundation closely.
- **Abandonment of a multi-phase settlement**: uniform within the last phase,
  or the phase end plus a floored exponential delay (`exp`).
- **Peaks**: uniform, or a floored truncated normal centred on its span with
  sd `peak_sd_perc` × span (`norm`). The first peak spans the foundation to the
  first phase end, the last spans the last phase start to the abandonment, and
  the others span their phase. A single-phase normal peak spans the foundation
  to the phase end, so it can fall after the abandonment.
- **Sizes**: 0.1 at foundation, the phase size at each peak, the last phase
  size at the end, or 0 when the end falls before `terminus`.
- **Growth**: R's `growth_rate()` orders each draw by date, computes
  `cagr = (size / previous) ** (1 / max(1, years)) - 1` and the linear rate,
  sets infinite values to missing, and writes rates and phase bounds back by
  position. When a single-phase normal peak follows the abandonment, the two
  rows therefore carry each other's rate, exactly as in R.
- **Site summaries** group by settlement and phase bounds, so the foundation
  joins the first phase and the abandonment the last. `median`, `mean` and
  `med_lin` are missing for a group containing a missing rate, as R's defaults
  leave them, which blanks the first phase. Rows are ordered as dplyr orders
  them, latest phase first.
- **Period summaries** skip missing values, and `size` is the sum over every
  row and draw, as R's `area`.
- **Persistence** uses only rows from phases starting before `terminus`. A
  settlement occupied in a phase ending at `terminus` gets the full
  `available_pers`, otherwise the span between its earliest and latest
  remaining date; a persistence of zero becomes one year after `norm_pers` is
  computed.
- **Plots** follow R: `plot_site_growth()` leaves the earliest phase out of the
  bands and y range, and so cannot draw a single-phase site;
  `plot_trajectories()` defaults to `growth_cagr`, `seed=1234` and R's y range
  (the data range widened by a tenth of each limit); both pad axes by 4 percent.

Summary columns keep R's names and units: `med_year` in BC/AD, and
`q025_year`, `q975_year` and the persistence `start` and `end` in years BP.
`year_median` is added in astronomical years for plotting.

## What necessarily differs

- **Random numbers.** Draws come from NumPy, not R's generator, so individual
  draws and seeds differ while the distributions are the same.
- **Truncated draws.** R redraws until a value falls inside its bounds and stops
  after 1000 failures. Here the truncated distribution is drawn directly, which
  gives the same distribution of accepted values but never fails. A
  zero-length exponential span, on which R stops with an error, gives a delay
  of zero.
- **Years in the trajectory frame** are astronomical, like the rest of the
  package. Convert with `mc.to_bp()`.
- **Site identifiers** are strings. Summaries order them numerically when all
  are numbers and as text otherwise, as R does for integer or character IDs.

## Running on the example datasets

[`scripts/run_sampling.py`](../scripts/run_sampling.py) samples Yautepec and
Viabundus with 500 iterations and writes figures and summary tables to
`figures/<dataset>/sampling-python/`, with the tables in the column layout of
the R output (`start` and `end` in years BP). For Yautepec it repeats
`scripts/r/run-sampling-plots.R`: the same seven sites, a second run of sites
160 and 200 with exponential dates and normal peaks, and the same file names,
so `figures/yautepec/sampling/` and `sampling-python/` compare file by file:

```bash
python scripts/run_sampling.py
python scripts/run_sampling.py --dataset viabundus
```

Viabundus is snapshot data, which `sample_trajectories()` rejects. The runner
reads each census as the period inferred around it (1300 as 1250-1350, up to
1650 as 1625-1675), drops zero-population rows and sorts each city's rows
oldest first, since phases are read in row order. The terminus is 1625 CE for
sampling and 1675 CE for persistence, the start and end of the last period.

Read the Viabundus results with three limits in mind. Census years are known
dates, so sampling a date within each period only spreads the estimates.
Cities seen at the first census existed before it, so the growth from the
seed size (`initial_size`, here 0.1 inhabitants) into the first period is an
artefact. It dominates the first period in `plot_period_growth()`, while
`plot_site_growth()` leaves it out as R does. And 377 of 400 cities are
present at every census, so persistence mostly measures the record's length.

## Verification against R

With 4000 draws for 13 Yautepec sites (the seven multi-phase focus sites and
six single-phase ones) under three method combinations, Python and R agree on
every draw position. Two-sample KS tests on the 291 date distributions show no
differences beyond chance, the share of draws abandoned to zero and of missing
growth rates match, and `summary_site()`, `summary_period()` and
`persistence()` agree within Monte Carlo noise, including which medians are
missing.
