"""Right-censored lifespans after chronological sampling.

:func:`modelcity.sampling.persistence_samples` stops at the end of the record.
A settlement still occupied then has ``norm_pers`` of 1: the draw says how long
it had already lasted, not how long it would have lasted.  The functions here
treat that as a right-censored survival time, fit the parametric families used
by ``flexsurv``, and draw a total lifespan for the censored rows from the
generalised gamma.

This follows ``scripts/r/longevity_functions.R``.  The hazard shapes are those
compared by Scheffer et al. (PNAS 2023, ``10.1073/pnas.2218834120``).  Fits use
SciPy's optimiser rather than ``flexsurvreg``, so coefficients are close on a
shared sample and are not identical.  Imputed draws use NumPy's generator.
"""

from __future__ import annotations

import warnings
from typing import Dict, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from scipy import optimize, stats

from .sampling import _r_order

__all__ = [
    "SURVIVAL_DISTRIBUTIONS",
    "SurvivalFit",
    "mark_censoring",
    "fit_survival",
    "compare_survival",
    "fit_persistence",
    "impute_lifespan",
    "summarize_persistence",
]

#: Names accepted by ``flexsurvreg``.  ``llogis`` is the log-logistic.
SURVIVAL_DISTRIBUTIONS = (
    "exponential",
    "weibull",
    "lognormal",
    "gompertz",
    "gamma",
    "llogis",
    "gengamma",
)

_ALIASES = {
    "loglogistic": "llogis",
    "log-logistic": "llogis",
    "generalized gamma": "gengamma",
    "generalised gamma": "gengamma",
}

_PARAMETERS = {
    "exponential": ("rate",),
    "weibull": ("shape", "scale"),
    "lognormal": ("meanlog", "sdlog"),
    "gompertz": ("shape", "rate"),
    "gamma": ("shape", "rate"),
    "llogis": ("shape", "scale"),
    "gengamma": ("mu", "sigma", "Q"),
}

_SUMMARY_QUANTILES = (
    ("q025", 0.025),
    ("q975", 0.975),
    ("q25", 0.25),
    ("q75", 0.75),
)


class SurvivalFit:
    """One converged parametric fit on the scale reported by ``flexsurv``.

    ``parameters`` uses ``flexsurv``'s names (``rate``, ``shape`` and
    ``scale``, ``meanlog`` and ``sdlog``, or ``mu``, ``sigma`` and ``Q``).
    ``aic`` is ``-2 * log_likelihood + 2 * n_params``.
    """

    def __init__(self, distribution: str, parameters: Mapping[str, float], log_likelihood: float):
        self.distribution = distribution
        self.parameters = {name: float(parameters[name]) for name in _PARAMETERS[distribution]}
        self.log_likelihood = float(log_likelihood)
        self.n_params = len(self.parameters)
        self.aic = -2.0 * self.log_likelihood + 2.0 * self.n_params

    def __repr__(self) -> str:
        terms = ", ".join("{}={:.4g}".format(name, value) for name, value in self.parameters.items())
        return "SurvivalFit({}: {}, AIC={:.4g})".format(self.distribution, terms, self.aic)


def _prentice_gengamma(mu: float, sigma: float, q: float):
    """Prentice generalised gamma as a SciPy distribution.

    ``flexsurv``'s ``(mu, sigma, Q)`` converts to Stacy's ``(shape, k, scale)``
    by ``shape = 1/Q^2``, ``k = Q/sigma`` and
    ``scale = exp(mu + sigma * log(Q^2) / Q)``.  ``Q = 0`` is the lognormal
    limit, with ``meanlog = mu`` and ``sdlog = sigma``.
    """
    if abs(q) < 1e-8:
        return stats.lognorm(sigma, scale=np.exp(mu))
    shape = 1.0 / q**2
    power = q / sigma
    scale = np.exp(mu + (sigma / q) * np.log(q * q))
    return stats.gengamma(shape, power, scale=scale)


def _require_lifespans(longevity) -> np.ndarray:
    times = np.asarray(longevity, dtype=float).ravel()
    if times.size == 0 or not np.isfinite(times).all() or np.any(times <= 0):
        raise ValueError("Lifespans must be strictly positive and finite.")
    return times


def _require_event(event, n_times: int) -> np.ndarray:
    status = np.asarray(event, dtype=float).ravel()
    if status.shape != (n_times,):
        raise ValueError("event must have one value for each lifespan.")
    if not np.isfinite(status).all() or not np.isin(status, (0.0, 1.0)).all():
        raise ValueError("event must be 0 (still occupied) or 1 (abandonment observed).")
    return status


def _resolve_distributions(distributions: Union[str, Sequence[str]]) -> List[str]:
    if distributions == "all":
        return list(SURVIVAL_DISTRIBUTIONS)
    requested = [distributions] if isinstance(distributions, str) else list(distributions)
    if not requested:
        raise ValueError("distributions must be 'all' or at least one distribution name.")
    resolved = []
    for name in requested:
        canonical = _ALIASES.get(str(name).lower(), str(name).lower())
        if canonical not in _PARAMETERS:
            raise ValueError(
                "unknown distribution {!r}. Choose 'all' or any of: {}.".format(
                    name, ", ".join(SURVIVAL_DISTRIBUTIONS)
                )
            )
        if canonical not in resolved:
            resolved.append(canonical)
    return resolved


def _censored_nll(log_density: np.ndarray, log_survival: np.ndarray, event: np.ndarray) -> float:
    log_likelihood = np.where(event == 1.0, log_density, log_survival)
    if not np.isfinite(log_likelihood).all():
        return np.inf
    return -float(log_likelihood.sum())


def _minimize(objective, starts: Sequence[Sequence[float]]):
    best = None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for start in starts:
            try:
                found = optimize.minimize(objective, np.asarray(start, dtype=float), method="L-BFGS-B")
            except (ValueError, FloatingPointError, RuntimeError):
                continue
            if not np.isfinite(found.fun):
                continue
            if best is None or found.fun < best.fun:
                best = found
    return best


def _location_scale(times: np.ndarray) -> Tuple[float, float]:
    logged = np.log(times)
    centre = float(np.median(logged))
    spread = float(np.std(logged, ddof=1)) if times.size > 1 else 1.0
    if not np.isfinite(spread) or spread < 0.05:
        spread = 0.5
    return centre, spread


def _fit_exponential(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    def objective(point):
        rate = np.exp(point[0])
        return _censored_nll(np.log(rate) - rate * times, -rate * times, event)

    found = _minimize(objective, [[np.log(1.0 / times.mean())]])
    if found is None:
        return None
    return SurvivalFit("exponential", {"rate": np.exp(found.x[0])}, -found.fun)


def _fit_weibull(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    def objective(point):
        shape, scale = np.exp(point)
        distribution = stats.weibull_min(shape, scale=scale)
        return _censored_nll(distribution.logpdf(times), distribution.logsf(times), event)

    found = _minimize(objective, [[0.0, np.log(times.mean())]])
    if found is None:
        return None
    shape, scale = np.exp(found.x)
    return SurvivalFit("weibull", {"shape": shape, "scale": scale}, -found.fun)


def _fit_lognormal(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    centre, spread = _location_scale(times)

    def objective(point):
        meanlog, sdlog = point[0], np.exp(point[1])
        distribution = stats.lognorm(sdlog, scale=np.exp(meanlog))
        return _censored_nll(distribution.logpdf(times), distribution.logsf(times), event)

    found = _minimize(objective, [[centre, np.log(spread)]])
    if found is None:
        return None
    return SurvivalFit(
        "lognormal",
        {"meanlog": found.x[0], "sdlog": np.exp(found.x[1])},
        -found.fun,
    )


def _fit_gompertz(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    def objective(point):
        shape, rate = point[0], np.exp(point[1])
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            if abs(shape) < 1e-8:
                cumulative = rate * times
                log_hazard = np.full(times.shape, np.log(rate))
            else:
                cumulative = rate / shape * (np.exp(shape * times) - 1.0)
                log_hazard = np.log(rate) + shape * times
            log_survival = -cumulative
            log_density = log_hazard + log_survival
        return _censored_nll(log_density, log_survival, event)

    rate0 = 1.0 / times.mean()
    found = _minimize(objective, [[0.0, np.log(rate0)], [-0.01, np.log(rate0)]])
    if found is None:
        return None
    return SurvivalFit("gompertz", {"shape": found.x[0], "rate": np.exp(found.x[1])}, -found.fun)


def _fit_gamma(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    def objective(point):
        shape, rate = np.exp(point)
        distribution = stats.gamma(shape, scale=1.0 / rate)
        return _censored_nll(distribution.logpdf(times), distribution.logsf(times), event)

    found = _minimize(objective, [[0.0, np.log(1.0 / times.mean())]])
    if found is None:
        return None
    shape, rate = np.exp(found.x)
    return SurvivalFit("gamma", {"shape": shape, "rate": rate}, -found.fun)


def _fit_loglogistic(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    def objective(point):
        shape, scale = np.exp(point)
        distribution = stats.fisk(shape, scale=scale)
        return _censored_nll(distribution.logpdf(times), distribution.logsf(times), event)

    found = _minimize(objective, [[0.0, np.log(np.median(times))]])
    if found is None:
        return None
    shape, scale = np.exp(found.x)
    return SurvivalFit("llogis", {"shape": shape, "scale": scale}, -found.fun)


def _fit_gengamma(times: np.ndarray, event: np.ndarray) -> Optional[SurvivalFit]:
    centre, spread = _location_scale(times)

    def objective(point):
        mu, sigma, q = point[0], np.exp(point[1]), point[2]
        if abs(q) < 1e-4:
            return np.inf
        distribution = _prentice_gengamma(mu, sigma, q)
        return _censored_nll(distribution.logpdf(times), distribution.logsf(times), event)

    starts = [[centre, np.log(spread), sign] for sign in (1.0, -1.0)]
    found = _minimize(objective, starts)
    if found is None:
        return None
    mu, sigma, q = found.x[0], np.exp(found.x[1]), found.x[2]
    return SurvivalFit("gengamma", {"mu": mu, "sigma": sigma, "Q": q}, -found.fun)


_FITTERS = {
    "exponential": _fit_exponential,
    "weibull": _fit_weibull,
    "lognormal": _fit_lognormal,
    "gompertz": _fit_gompertz,
    "gamma": _fit_gamma,
    "llogis": _fit_loglogistic,
    "gengamma": _fit_gengamma,
}


def mark_censoring(samples: pd.DataFrame) -> pd.DataFrame:
    """Add the event indicator used by a survival fit.

    ``event`` is 0 when ``norm_pers`` is 1, the settlement still being occupied
    at the terminus, and 1 when the abandonment was observed.  This is the
    ``survival`` column in ``longevity.R``: 0 is censored and 1 is the event,
    which is the coding ``Surv`` expects.  A missing ``norm_pers`` stays missing.
    """
    if "norm_pers" not in samples.columns:
        raise ValueError("samples must contain norm_pers, as returned by persistence_samples().")
    marked = samples.copy()
    ratio = marked["norm_pers"].to_numpy(dtype=float)
    event = np.full(ratio.shape, np.nan)
    known = np.isfinite(ratio)
    event[known] = np.where(ratio[known] == 1.0, 0.0, 1.0)
    marked["event"] = event
    return marked


def fit_survival(
    longevity,
    event,
    distributions: Union[str, Sequence[str]] = "all",
) -> Dict[str, SurvivalFit]:
    """Fit right-censored parametric models, as ``fit_functions()``.

    ``longevity`` is strictly positive years.  ``event`` is 1 when the
    abandonment was observed and 0 when the settlement was still occupied.
    ``distributions`` is ``"all"`` or any of :data:`SURVIVAL_DISTRIBUTIONS`
    (``loglogistic`` is accepted for ``llogis``).  A distribution that does
    not converge is omitted, as ``flexsurvreg`` failures are in R.
    """
    times = _require_lifespans(longevity)
    status = _require_event(event, times.size)
    fitted = {}
    for name in _resolve_distributions(distributions):
        fit = _FITTERS[name](times, status)
        if fit is not None and np.isfinite(fit.log_likelihood):
            fitted[name] = fit
    return fitted


def compare_survival(
    longevity,
    event,
    distributions: Union[str, Sequence[str]] = "all",
) -> pd.DataFrame:
    """Rank fitted models by AIC, as ``compare_functions()``.

    Columns follow R: ``Model``, ``AIC``, ``N_params``, ``Delta_AIC``,
    ``Weight`` and ``Rank``.  Weights are Akaike weights from the fitted set.
    """
    fitted = fit_survival(longevity, event, distributions)
    if not fitted:
        raise ValueError("no distribution converged.")
    comparison = pd.DataFrame(
        {
            "Model": list(fitted),
            "AIC": [fit.aic for fit in fitted.values()],
            "N_params": [fit.n_params for fit in fitted.values()],
        }
    )
    comparison = comparison.sort_values(["AIC", "Model"], kind="mergesort").reset_index(drop=True)
    delta = comparison["AIC"] - comparison["AIC"].min()
    weight = np.exp(-0.5 * delta)
    comparison["Delta_AIC"] = delta
    comparison["Weight"] = weight / weight.sum()
    comparison["Rank"] = np.arange(1, len(comparison) + 1)
    return comparison


def fit_persistence(samples: pd.DataFrame, distribution: str = "gengamma") -> pd.DataFrame:
    """Fit one distribution to each Monte Carlo iteration.

    ``samples`` is a :func:`modelcity.sampling.persistence_samples` table.
    ``event`` is added from ``norm_pers`` when it is not already present.
    The result has one row per iteration that converged, with that
    distribution's parameters plus ``loglik`` and ``aic``.  Averaging the
    parameter columns is the step ``longevity.R`` takes before imputation.
    """
    names = _resolve_distributions(distribution)
    if len(names) != 1:
        raise ValueError("fit_persistence fits one distribution; pass fit_survival the others.")
    name = names[0]
    frame = samples if "event" in samples.columns else mark_censoring(samples)
    if "iteration" not in frame.columns or "persistence" not in frame.columns:
        raise ValueError("samples must contain iteration and persistence.")

    rows = []
    for iteration, group in frame.groupby("iteration", sort=True):
        usable = group["persistence"].notna() & group["event"].notna()
        subset = group.loc[usable]
        if subset.empty:
            continue
        fitted = fit_survival(subset["persistence"], subset["event"], distributions=name)
        fit = fitted.get(name)
        if fit is None:
            continue
        rows.append({"iteration": iteration, **fit.parameters, "loglik": fit.log_likelihood, "aic": fit.aic})
    columns = ["iteration", *_PARAMETERS[name], "loglik", "aic"]
    return pd.DataFrame(rows, columns=columns)


def impute_lifespan(longevity, mu: float, sigma: float, q: float, seed: int = 1234) -> np.ndarray:
    """Draw a total lifespan from the generalised gamma, conditional on surviving so far.

    For each observed duration ``t`` the draw is ``floor(F^{-1}(F(t) + U(1-F(t))))``,
    the total time rather than the remaining years.  ``mu``, ``sigma`` and ``q``
    are the Prentice parameters and are the ones used for both the cdf and the
    quantile.  A non-finite quantile becomes missing.  ``seed`` feeds NumPy's
    generator, so the same seed repeats and does not match an R ``set.seed``.
    """
    times = _require_lifespans(longevity)
    if not np.isfinite(sigma) or sigma <= 0 or not np.isfinite(mu) or not np.isfinite(q):
        raise ValueError("mu and Q must be finite, and sigma must be positive.")
    distribution = _prentice_gengamma(float(mu), float(sigma), float(q))
    uniform = np.random.default_rng(seed).random(times.size)
    cdf = distribution.cdf(times)
    conditional = cdf + uniform * (1.0 - cdf)
    with np.errstate(over="ignore", invalid="ignore"):
        drawn = distribution.ppf(conditional)
        total = np.floor(drawn)
    total[~np.isfinite(total)] = np.nan
    below = np.isfinite(total) & (total < times)
    total[below] = np.ceil(times[below])
    return total


def summarize_persistence(samples: pd.DataFrame) -> pd.DataFrame:
    """Spread of an existing per-draw persistence table, as ``long_persistence()``.

    Unlike :func:`modelcity.sampling.persistence_summary`, this does not
    recompute spans from trajectory dates, so it can summarise a table whose
    censored lifespans have been replaced.  The ``_n`` columns average the
    ``norm_pers`` already stored on the rows.  ``long_persistence()`` does the
    same, so a share left at its pre-imputation value stays there: censored
    rows remain 1, and the mean share stays within 0 and 1.  Columns match
    ``persistence_summary``.
    """
    required = ("id", "persistence", "norm_pers", "start", "end")
    missing = [column for column in required if column not in samples.columns]
    if missing:
        raise ValueError("samples is missing {}.".format(", ".join(missing)))

    grouped = samples.groupby("id", sort=False)
    result = pd.DataFrame(
        {"mean": grouped["persistence"].mean(), "median": grouped["persistence"].median()}
    )
    for name, probability in _SUMMARY_QUANTILES:
        result[name] = grouped["persistence"].quantile(probability)
    result["mean_n"] = grouped["norm_pers"].mean()
    result["median_n"] = grouped["norm_pers"].median()
    for name, probability in _SUMMARY_QUANTILES:
        result[name + "_n"] = grouped["norm_pers"].quantile(probability)
    result["start"] = grouped["start"].median()
    result["end"] = grouped["end"].median()
    return result.reset_index().sort_values("id", key=_r_order).reset_index(drop=True)
