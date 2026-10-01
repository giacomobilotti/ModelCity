"""Monte Carlo sampling of when settlements were founded, peaked and ended.

A period only bounds when something happened: a site recorded in a 400-year
phase could have been founded, reached its recorded size, or been abandoned at
any point inside it.  The functions here draw those dates repeatedly within the
phase bounds, pair each date with the recorded size, and carry every draw
through to growth rates and persistence.

This reproduces ``scripts/r/sampling.R`` and ``scripts/r/helpers_sampling.R``:
the same rules, whole-year flooring, ``terminus`` handling and summary
statistics.  Draws come from NumPy rather than R's generator, so individual
draws differ while their distributions are the same.  The rules are applied in
years BP, as in R, and years are returned as astronomical calendar years, like
:attr:`Settlements.data`.  Summary columns that share a name with an R output
column also share its units.

Every draw of one settlement is a sequence of *steps*:

* step 0, ``foundation``, in the first phase, at ``initial_size``;
* one ``peak`` per phase, at that phase's recorded size;
* a final ``abandonment``, in the last phase, at the last recorded size, or
  at zero when it falls before ``terminus``.
"""

from __future__ import annotations

from typing import Optional, Union

import numpy as np
import pandas as pd
from scipy.special import ndtr, ndtri

from .schema import SchemaError
from .settlements import Settlements
from .timescale import bp_to_ce, format_period, to_astronomical, to_bp

__all__ = [
    "START_METHODS",
    "END_METHODS",
    "PEAK_METHODS",
    "sample_trajectories",
    "site_summary",
    "period_summary",
    "persistence_samples",
    "persistence_summary",
]

START_METHODS = ("unif", "exp")
END_METHODS = ("unif", "exp")
PEAK_METHODS = ("unif", "norm")

_ALIASES = {"uniform": "unif", "exponential": "exp", "normal": "norm"}
_EVENTS = ("foundation", "peak", "abandonment")

Seed = Union[None, int, np.random.Generator]


def _method(value: str, allowed, role: str) -> str:
    resolved = _ALIASES.get(value, value)
    if resolved not in allowed:
        raise ValueError(
            "{} must be one of {}, got {!r}.".format(
                role, " or ".join(repr(option) for option in allowed), value
            )
        )
    return resolved


def _terminus_bp(terminus: float, time_scale: str) -> float:
    return float(terminus) if time_scale == "BP" else float(to_bp(terminus))


# -- samplers, all in years BP ----------------------------------------------

def _floor_uniform(rng: np.random.Generator, low, high, n: int) -> np.ndarray:
    """``floor(runif(1, min = low, max = high))``."""
    return np.floor(low + (high - low) * rng.random(n))


def _exp_delay(rng: np.random.Generator, delta, perc: float, n: int) -> np.ndarray:
    """``exp_sampling()``: ``floor(rexp(1, perc / delta))``, redrawn until below ``delta``.

    Drawn by inverting the truncated CDF, which gives the distribution of the
    accepted draws without rejecting any.
    """
    delta = np.broadcast_to(np.asarray(delta, dtype=float), (n,))
    # floor(x) < delta exactly when x < ceil(delta).
    bound = np.ceil(delta)
    with np.errstate(divide="ignore", invalid="ignore"):
        rate = perc / delta
        offset = -np.log1p(-rng.random(n) * -np.expm1(-rate * bound)) / rate
    delay = np.minimum(np.floor(offset), bound - 1)
    # R cannot sample a zero-length span at all and stops; no delay is the
    # only value that keeps the date inside it.
    return np.where(delta > 0, delay, 0.0)


def _norm_peak(rng: np.random.Generator, start, end, peak_sd_perc: float, n: int) -> np.ndarray:
    """``norm_peak()``: ``floor(rnorm())`` around the midpoint, redrawn until within ``[end, start]``."""
    start = np.broadcast_to(np.asarray(start, dtype=float), (n,))
    end = np.broadcast_to(np.asarray(end, dtype=float), (n,))
    phase = start - end
    mean = start - phase / 2.0
    if peak_sd_perc == 0:
        return np.floor(mean)

    sd = phase * peak_sd_perc
    # floor(y) lands in [end, start] exactly when y is in [ceil(end), floor(start) + 1).
    low = np.ceil(end)
    high = np.floor(start) + 1.0
    with np.errstate(divide="ignore", invalid="ignore"):
        lower_p = ndtr((low - mean) / sd)
        upper_p = ndtr((high - mean) / sd)
        draw = mean + sd * ndtri(lower_p + rng.random(n) * (upper_p - lower_p))
    peak = np.clip(np.floor(draw), low, high - 1.0)
    return np.where(sd > 0, peak, np.floor(mean))


# -- sampling ---------------------------------------------------------------

def sample_trajectories(
    settlements: Settlements,
    iterations: int = 1000,
    *,
    start_method: str = "unif",
    end_method: str = "unif",
    peak_method: str = "unif",
    perc: float = 5.0,
    peak_sd_perc: float = 0.1,
    terminus: Optional[float] = None,
    initial_size: float = 0.1,
    min_interval: float = 1.0,
    seed: Seed = None,
) -> pd.DataFrame:
    """Draw ``iterations`` dated size trajectories for every settlement.

    This is ``run_mc()``: ``sample_dates()`` followed by ``growth_rate()``,
    repeated.  Each row of a settlement is one phase, taken in input order.

    Parameters
    ----------
    settlements:
        Interval records, at most one row per settlement and phase start.
    iterations:
        Number of draws per settlement.
    start_method:
        ``"unif"`` places the foundation anywhere in the first phase;
        ``"exp"`` favours its beginning.
    end_method:
        ``"unif"`` places the abandonment anywhere in the last phase.
        ``"exp"`` favours the end of the last phase for a multi-phase
        settlement, and the foundation date for a single-phase one.
    peak_method:
        ``"unif"`` or ``"norm"``.  A normal peak is centred on its span, with
        a standard deviation of ``peak_sd_perc`` times the span.  For a
        single-phase settlement that span runs from the foundation to the
        phase end, so the peak can fall after the sampled abandonment.
    perc:
        Shape of the exponential methods: before truncation the mean offset
        is the span divided by ``perc``.
    peak_sd_perc:
        Standard deviation of a normal peak, as a fraction of its span.
    terminus:
        A settlement whose sampled abandonment falls before this year ends at
        size zero.  Given in the dataset's own time scale (years BP for a BP
        dataset).  Defaults to 1 BP, as in R, so every settlement that ends
        before 1949 CE ends at zero.
    initial_size:
        Size at foundation.
    min_interval:
        Shortest span, in years, used as the denominator of a growth rate.
    seed:
        Seed or generator, for reproducible draws.

    Returns
    -------
    pandas.DataFrame
        One row per settlement, draw and step, with ``id``, ``iteration``,
        ``step``, ``event``, ``phase``, phase bounds ``t_start`` and
        ``t_end``, the sampled ``year``, ``size``, and the growth from the
        previous date, ``growth_cagr`` and ``growth_linear``.

    Notes
    -----
    As in R's ``growth_rate()``, growth is computed after ordering each draw
    by date and written back by position, together with the phase bounds.
    When a single-phase normal peak falls after the abandonment, the last two
    rows of that draw therefore carry each other's growth rate.
    """
    start_method = _method(start_method, START_METHODS, "start_method")
    end_method = _method(end_method, END_METHODS, "end_method")
    peak_method = _method(peak_method, PEAK_METHODS, "peak_method")
    if iterations < 1:
        raise ValueError("iterations must be at least 1, got {}.".format(iterations))
    if perc <= 0:
        raise ValueError("perc must be > 0, got {}.".format(perc))
    if peak_sd_perc < 0:
        raise ValueError("peak_sd_perc must be > 0, got {}.".format(peak_sd_perc))

    if settlements.is_snapshot:
        raise SchemaError(
            "sampling needs interval records with a start and an end per phase; "
            "snapshot observations have no phase to sample dates from."
        )

    data = settlements.data
    if data.empty:
        raise ValueError("no phases to sample.")
    empty = data["t_start"] >= data["t_end"]
    if empty.any():
        raise SchemaError(
            "start must be older than end for every phase.\n"
            "  sites affected: {}".format(", ".join(map(str, data.loc[empty, "id"].unique()[:10])))
        )
    repeated = data.duplicated(["id", "t_start"], keep=False)
    if repeated.any():
        raise SchemaError(
            "duplicates found: settlements with more than one row for the same phase start.\n"
            "  sites affected: {}".format(", ".join(map(str, data.loc[repeated, "id"].unique()[:10])))
        )

    terminus_bp = 1.0 if terminus is None else _terminus_bp(terminus, settlements.time.time_scale)
    rng = seed if isinstance(seed, np.random.Generator) else np.random.default_rng(seed)
    n = int(iterations)

    site_ids = []
    blocks = {key: [] for key in ("site", "step", "event", "phase", "start", "end", "year", "size")}

    for index, (site_id, phases) in enumerate(data.groupby("id", sort=False)):
        starts = to_bp(phases["t_start"].to_numpy(dtype=float))
        ends = to_bp(phases["t_end"].to_numpy(dtype=float))
        sizes = phases["size"].to_numpy(dtype=float)
        k = len(starts)

        if start_method == "unif":
            first = _floor_uniform(rng, ends[0], starts[0], n)
        else:
            first = starts[0] - _exp_delay(rng, starts[0] - ends[0], perc, n)

        if k == 1:
            if end_method == "unif":
                last = _floor_uniform(rng, ends[0], first, n)
            else:
                last = first - _exp_delay(rng, first - ends[0], perc, n)
            if peak_method == "unif":
                peaks = [_floor_uniform(rng, last, first, n)]
            else:
                peaks = [_norm_peak(rng, first, ends[0], peak_sd_perc, n)]
        else:
            if end_method == "unif":
                last = _floor_uniform(rng, ends[-1], starts[-1], n)
            else:
                last = ends[-1] + _exp_delay(rng, starts[-1] - ends[-1], perc, n)
            if peak_method == "unif":
                peaks = [_floor_uniform(rng, ends[0], first, n)]
                peaks += [_floor_uniform(rng, ends[i], starts[i], n) for i in range(1, k - 1)]
                peaks.append(_floor_uniform(rng, last, starts[-1], n))
            else:
                peaks = [_norm_peak(rng, first, ends[0], peak_sd_perc, n)]
                peaks += [_norm_peak(rng, starts[i], ends[i], peak_sd_perc, n) for i in range(1, k - 1)]
                peaks.append(_norm_peak(rng, starts[-1], last, peak_sd_perc, n))

        steps = k + 2
        years = np.column_stack([first] + peaks + [last])
        size_of_step = np.tile(np.concatenate([[initial_size], sizes, [sizes[-1]]]), (n, 1))
        size_of_step[:, -1] = np.where(last > terminus_bp, 0.0, sizes[-1])
        phase_of_step = np.concatenate([[0], np.arange(k), [k - 1]])

        site_ids.append(site_id)
        blocks["site"].append(np.full(n * steps, index))
        blocks["step"].append(np.tile(np.arange(steps), n))
        blocks["event"].append(np.tile([0] + [1] * k + [2], n))
        blocks["phase"].append(np.tile(phase_of_step, n))
        blocks["start"].append(np.tile(starts[phase_of_step], n))
        blocks["end"].append(np.tile(ends[phase_of_step], n))
        blocks["year"].append(years.ravel())
        blocks["size"].append(size_of_step.ravel())

    joined = {key: np.concatenate(values) for key, values in blocks.items()}
    steps_per_site = np.bincount(joined["site"]) // n
    frame = pd.DataFrame(
        {
            "id": pd.Categorical.from_codes(joined["site"], categories=site_ids),
            "iteration": np.concatenate(
                [np.repeat(np.arange(1, n + 1), steps) for steps in steps_per_site]
            ),
            "step": joined["step"],
            "event": pd.Categorical.from_codes(joined["event"], categories=list(_EVENTS)),
            "phase": joined["phase"],
            "t_start": to_astronomical(joined["start"]),
            "t_end": to_astronomical(joined["end"]),
            "year": to_astronomical(joined["year"]),
            "size": joined["size"],
        }
    )
    _growth_rate(frame, min_interval)
    frame.attrs.update(
        time_scale=settlements.time.time_scale,
        size_type=settlements.size.size_type,
        size_unit=settlements.size.unit,
        size_label=settlements.size.label,
    )
    return frame


def _growth_rate(frame: pd.DataFrame, min_interval: float) -> None:
    """R's ``growth_rate()`` for rows grouped by draw and ordered by step.

    Each draw is ordered by date, rates are computed on that order, and the
    rates and phase bounds are written back by position, exactly as R assigns
    a sorted subset into the unsorted frame.
    """
    first = frame["step"].to_numpy() == 0
    draw = np.cumsum(first)
    order = np.lexsort((frame["year"].to_numpy(), draw))

    years = frame["year"].to_numpy(dtype=float)[order]
    sizes = frame["size"].to_numpy(dtype=float)[order]
    previous_year = np.roll(years, 1)
    previous_size = np.roll(sizes, 1)
    previous_year[first] = np.nan
    previous_size[first] = np.nan

    interval = np.maximum(years - previous_year, min_interval)
    with np.errstate(divide="ignore", invalid="ignore"):
        cagr = (sizes / previous_size) ** (1.0 / interval) - 1.0
        linear = (sizes - previous_size) / interval
    cagr[np.isinf(cagr)] = np.nan
    linear[np.isinf(linear)] = np.nan

    frame["t_start"] = frame["t_start"].to_numpy()[order]
    frame["t_end"] = frame["t_end"].to_numpy()[order]
    frame["growth_cagr"] = cagr
    frame["growth_linear"] = linear


# -- summaries --------------------------------------------------------------

def _median_or_na(grouped, column: str) -> pd.Series:
    """R's ``median()`` without ``na.rm``: missing as soon as one value is."""
    median = grouped[column].median()
    complete = grouped[column].count() == grouped[column].size()
    return median.where(complete)


def _mean_or_na(grouped, column: str) -> pd.Series:
    mean = grouped[column].mean()
    complete = grouped[column].count() == grouped[column].size()
    return mean.where(complete)


def _quantile(grouped, column: str, q: float) -> pd.Series:
    return grouped[column].quantile(q)


def _growth_statistics(grouped, keep_missing: bool) -> pd.DataFrame:
    if keep_missing:
        median, mean = _median_or_na(grouped, "growth_cagr"), _mean_or_na(grouped, "growth_cagr")
        median_linear = _median_or_na(grouped, "growth_linear")
    else:
        median, mean = grouped["growth_cagr"].median(), grouped["growth_cagr"].mean()
        median_linear = grouped["growth_linear"].median()

    q25 = _quantile(grouped, "growth_cagr", 0.25)
    q75 = _quantile(grouped, "growth_cagr", 0.75)
    bp = grouped["year_bp"]
    return pd.DataFrame(
        {
            "median": median,
            "mean": mean,
            "q025": _quantile(grouped, "growth_cagr", 0.025),
            "q975": _quantile(grouped, "growth_cagr", 0.975),
            "q25": q25,
            "q75": q75,
            "iqr": q75 - q25,
            "med_lin": median_linear,
            "q025_lin": _quantile(grouped, "growth_linear", 0.025),
            "q975_lin": _quantile(grouped, "growth_linear", 0.975),
            "med_year": grouped["year_bcad"].median(),
            "q025_year": bp.quantile(0.025),
            "q975_year": bp.quantile(0.975),
            "year_median": grouped["year"].median(),
        }
    )


def _with_r_years(trajectories: pd.DataFrame) -> pd.DataFrame:
    frame = trajectories[["id", "t_start", "t_end", "year", "size", "growth_cagr", "growth_linear"]].copy()
    frame["year_bp"] = to_bp(frame["year"])
    frame["year_bcad"] = bp_to_ce(frame["year_bp"])
    return frame


def site_summary(trajectories: pd.DataFrame) -> pd.DataFrame:
    """R's ``summary_site()``: growth across draws per settlement and phase.

    Rows group by ``id`` and phase bounds, so the foundation joins the first
    phase and the abandonment the last.  ``median``, ``mean`` and ``med_lin``
    are missing when any value in the group is, as R's defaults leave them;
    the foundation therefore blanks them for the first phase.  Rows are
    ordered as in R, latest phase first.

    ``med_year`` is in BC/AD and ``q025_year`` and ``q975_year`` in years BP,
    the units R reports.  ``year_median`` is the astronomical median, used for
    plotting, and ``size`` is the largest size in the group.
    """
    frame = _with_r_years(trajectories)
    grouped = frame.groupby(["id", "t_start", "t_end"], observed=True, sort=False)
    result = _growth_statistics(grouped, keep_missing=True)
    result["size"] = grouped["size"].max()
    result = result.reset_index()
    result["id"] = result["id"].astype(str)
    return result.sort_values(["id", "t_start"], ascending=[True, False], key=_r_order).reset_index(drop=True)


def _r_order(column: pd.Series) -> pd.Series:
    """Sort keys as dplyr orders groups: numbers numerically, text by code point."""
    if column.name == "id":
        numeric = pd.to_numeric(column, errors="coerce")
        return numeric if numeric.notna().all() else column
    return column


def period_summary(trajectories: pd.DataFrame) -> pd.DataFrame:
    """R's ``summary_period()``: growth across settlements and draws per phase.

    Missing values are skipped, as R does here.  ``size`` is the sum over
    every row and draw, as R's ``area`` is.  Year columns follow
    :func:`site_summary`, and rows are ordered latest phase first.
    """
    frame = _with_r_years(trajectories)
    grouped = frame.groupby(["t_start", "t_end"], sort=False)
    result = _growth_statistics(grouped, keep_missing=False)
    result["size"] = grouped["size"].sum()
    result = result.reset_index().sort_values("t_start", ascending=False).reset_index(drop=True)
    result.insert(2, "label", [format_period(row.t_start, row.t_end) for row in result.itertuples()])
    return result


def _resolve_scale(trajectories: pd.DataFrame, time_scale: Optional[str]) -> str:
    scale = time_scale or trajectories.attrs.get("time_scale")
    if scale not in ("BP", "CE"):
        raise ValueError(
            "pass time_scale='BP' or 'CE' to say how terminus is expressed; "
            "the trajectories no longer carry it."
        )
    return scale


def persistence_samples(
    trajectories: pd.DataFrame,
    terminus: float,
    time_scale: Optional[str] = None,
) -> pd.DataFrame:
    """R's ``mc_persistence()``: how long each settlement lasted in each draw.

    Only rows from phases starting before ``terminus`` are used.  A
    settlement occupied in a phase ending at ``terminus`` persists for the
    whole ``available_pers``, from its foundation to ``terminus``; any other
    persists from its earliest to its latest remaining date.  ``norm_pers``
    is the share of the available time, and a persistence of zero is raised
    to one year after that share is taken.

    ``terminus`` is in the dataset's time scale, which ``time_scale``
    overrides if the frame has lost its metadata.  ``start`` and ``end`` are
    in years BP, as in R.
    """
    terminus_bp = _terminus_bp(terminus, _resolve_scale(trajectories, time_scale))
    frame = pd.DataFrame(
        {
            "id": trajectories["id"].astype(str),
            "iteration": trajectories["iteration"],
            "start_bp": to_bp(trajectories["t_start"]),
            "end_bp": to_bp(trajectories["t_end"]),
            "year_bp": to_bp(trajectories["year"]),
        }
    )
    frame = frame[frame["start_bp"] > terminus_bp]
    grouped = frame.groupby(["id", "iteration"], sort=True)
    oldest = grouped["year_bp"].max()
    youngest = grouped["year_bp"].min()
    reaches = np.isclose(grouped["end_bp"].min(), terminus_bp)

    available = oldest - terminus_bp
    persistence = available.where(reaches, oldest - youngest)
    normalised = persistence / available
    result = pd.DataFrame(
        {
            "available_pers": available,
            "persistence": persistence,
            "norm_pers": normalised,
            "start": oldest,
            "end": youngest.where(normalised != 1, terminus_bp),
        }
    )
    result.loc[result["persistence"] == 0, "persistence"] = 1.0
    return result.reset_index()


def persistence_summary(
    trajectories: pd.DataFrame,
    terminus: float,
    time_scale: Optional[str] = None,
) -> pd.DataFrame:
    """R's ``persistence()``: the spread of persistence across draws per settlement.

    Columns follow R: ``mean`` to ``q75`` in years, the same with an ``_n``
    suffix for the normalised share, and the median ``start`` and ``end`` in
    years BP.
    """
    samples = persistence_samples(trajectories, terminus, time_scale)
    grouped = samples.groupby("id", sort=False)
    result = pd.DataFrame({"mean": grouped["persistence"].mean(), "median": grouped["persistence"].median()})
    for name, q in (("q025", 0.025), ("q975", 0.975), ("q25", 0.25), ("q75", 0.75)):
        result[name] = grouped["persistence"].quantile(q)
    result["mean_n"] = grouped["norm_pers"].mean()
    result["median_n"] = grouped["norm_pers"].median()
    for name, q in (("q025_n", 0.025), ("q975_n", 0.975), ("q25_n", 0.25), ("q75_n", 0.75)):
        result[name] = grouped["norm_pers"].quantile(q)
    result["start"] = grouped["start"].median()
    result["end"] = grouped["end"].median()
    return result.reset_index().sort_values("id", key=_r_order).reset_index(drop=True)
