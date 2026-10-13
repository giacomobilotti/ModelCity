# Longevity

Persistence stops at the end of the record. A settlement still occupied at the
terminus has lasted exactly until that date in the draw, which is not the same
as a measured abandonment. Longevity treats that draw as right-censored, fits
parametric survival models to the durations, and samples a total lifespan for
the censored rows from a generalised gamma.

This is one of the three ModelCity analyses. It follows
[chronological sampling](chronological-sampling.md); the system-level figures
are in [Metrics and plots](metrics-and-plots.md). The functions live in
`modelcity.longevity`, are exported from `modelcity`, and reproduce the
reusable part of
[`scripts/r/longevity_functions.R`](../scripts/r/longevity_functions.R). They
take the table from `persistence_samples()`; snapshot data never reaches them
because `sample_trajectories()` rejects it. The hazard families are those
compared by Scheffer et al., *The vulnerability of aging states* (PNAS 2023,
[10.1073/pnas.2218834120](https://doi.org/10.1073/pnas.2218834120)).

[`scripts/run_longevity.py`](../scripts/run_longevity.py) runs the analysis on
the example datasets and writes `figures/<dataset>/longevity-python/`. It
runs Yautepec unless `--dataset` names another survey:

```bash
python scripts/run_longevity.py
python scripts/run_longevity.py --dataset viabundus
```

Settings for those runs are under
[Running on the example datasets](#running-on-the-example-datasets).

```python
import modelcity as mc

trajectories = mc.sample_trajectories(settlements, 1000, terminus=510, seed=42)
samples = mc.mark_censoring(mc.persistence_samples(trajectories, terminus=431))

fits = mc.fit_persistence(samples, distribution="gengamma")
coefficients = fits[["mu", "sigma", "Q"]].mean()

# Durations before imputation. The bins match geom_histogram(bins = 40).
mc.save_figure(mc.plot_longevity(samples, coefficients=coefficients), "outputs/longevity.png")

# Every draw of a site that was censored in any iteration, as in longevity.R.
# norm_pers is left unchanged, so the share in the summary is pre-imputation.
still_alive = samples.loc[samples["event"] == 0, "id"].unique()
replace = samples["id"].isin(still_alive)
samples.loc[replace, "persistence"] = mc.impute_lifespan(
    samples.loc[replace, "persistence"],
    coefficients["mu"], coefficients["sigma"], coefficients["Q"],
    seed=1234,
)
summary = mc.summarize_persistence(samples)
```

`compare_survival(longevity, event)` ranks the seven families by AIC on one
sample. `fit_survival` returns the fits themselves.

## Correspondence with R

| R | Python |
| --- | --- |
| `survival` column (`0` censored, `1` event) | `mark_censoring()` adds `event` with the same coding |
| `fit_functions()` | `fit_survival()` |
| `compare_functions()` | `compare_survival()` |
| per-iteration `gengamma` fit, then `colMeans` | `fit_persistence()`, then the mean of the parameter columns |
| `sample_remaining()` | `impute_lifespan()` |
| `long_persistence()` | `summarize_persistence()` |
| histogram with `dgengamma` | `plot_longevity()` |

`event` follows `Surv`: 1 means the abandonment was observed, 0 means the
settlement was still occupied. In the R script the column is called `survival`.

## What is reproduced

- **Censoring.** `norm_pers == 1` is censored. Any other finite share is an
  observed end. A missing share stays missing and is left out of a
  per-iteration fit.
- **Families.** Exponential, Weibull, lognormal, Gompertz, gamma, log-logistic
  (`llogis`; `loglogistic` is accepted too), and generalised gamma. A family
  that does not converge is omitted. Names and parameter scales match
  `flexsurvreg`: `rate`, `shape` and `scale`, `meanlog` and `sdlog`, and
  Prentice `mu`, `sigma`, `Q`.
- **Comparison.** Models are ordered by AIC. `Delta_AIC`, Akaike `Weight` and
  `Rank` follow `compare_functions()`.
- **Imputation.** For a censored duration `t`, draw
  `floor(F^{-1}(F(t) + U (1 - F(t))))` from the generalised gamma. The value
  is the total lifespan, not the remaining years. `mu`, `sigma` and `Q` are
  used for both the cdf and the quantile. A non-finite quantile is missing.
- **Summary.** `summarize_persistence()` pools an existing per-draw table by
  settlement, as `long_persistence()` does. Columns match
  `persistence_summary()`: `mean` through `q75` in years, the same with an
  `_n` suffix for the stored `norm_pers`, and the median `start` and `end` in
  years BP. The share is the column already on the rows. Leave it unchanged
  when imputing and the censored rows stay at 1, which is what the R script
  plots.

`plot_longevity()` draws a density histogram of the durations in the table. The
forty bins use the same edges as `ggplot2::geom_histogram(bins = 40)`, so the
first edge sits half a bin below the shortest duration. Pass the mean Prentice
coefficients and it adds the generalised-gamma density. That curve is the
density of a latent lifespan. The bars also include draws that were still
occupied at the terminus, so the curve can sit above the bars near zero and
miss a pile of long censored durations.

## What necessarily differs

- **Random numbers.** Imputation uses NumPy's generator. The default `seed`
  is 1234, as in R, but the same seed does not reproduce R's draws.
- **Optimiser.** Fits maximise the censored log-likelihood with SciPy, not
  `flexsurvreg`. On a shared sample the coefficients and AIC are close. They
  are not required to match digit for digit.
- **Arguments.** `impute_lifespan()` uses the `mu`, `sigma` and `Q` it is
  given. `sample_remaining()` evaluates the quantile from a global `best_c`.
- **Which rows are replaced.** `impute_lifespan()` changes only the durations
  it is given. The example runner, like `longevity.R`, passes every draw of a
  site that was censored in any iteration, including draws in which that site
  had already ended.

## What stays in the R script

`longevity.R` continues into a Yautepec analysis that the library does not
carry: Large / Failed Large / Small labels (10 ha, and a phase starting at
510 BP or later), soil zones, distance to rivers, distance to a site of at
least 10 ha (`dist_func()`), stepwise regressions, and a discrete-time hazard
for surviving into the next phase. Those steps need `yautepec.gpkg`, which is
not in this repository. The size labels are applied in
[`scripts/run_longevity.py`](../scripts/run_longevity.py) only so the example
figures can be compared with the R run.

## Running on the example datasets

[`scripts/run_longevity.py`](../scripts/run_longevity.py) samples every
settlement for 1000 iterations and writes figures and tables to
`figures/<dataset>/longevity-python/`:

```bash
python scripts/run_longevity.py
python scripts/run_longevity.py --dataset viabundus
```

Yautepec uses a sampling terminus of 510 BP and measures persistence to 431 BP.
Viabundus is snapshot data, so the runner first reads each census as the period
inferred around it, as [Chronological sampling](chronological-sampling.md)
describes, then samples to 1625 CE and measures persistence to 1675 CE. The
size-class figures are written for Yautepec only.

[`scripts/r/run-longevity-plots.R`](../scripts/r/run-longevity-plots.R) writes
the R counterparts to `figures/yautepec/longevity/`. Coefficients, the
normalised share (the pre-imputation `norm_pers`), and the histogram edges
follow the R run. Persistence histograms use `geom_histogram(bins = 40)`.
The draw from the mean generalised gamma uses base R `hist(breaks = 40)`,
which places breaks with `pretty()`. Individual imputed lifespans still
differ, because the two generators do.

## Verification against R

The Prentice generalised gamma agrees with `flexsurv`'s `dgengamma`,
`pgengamma` and `qgengamma` for positive and negative `Q`. On one fixed
censored sample, all seven families agree with `flexsurvreg` in their
parameter estimates, log-likelihood and AIC, within the tolerance of the two
optimisers.
