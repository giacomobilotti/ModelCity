"""Analytical summaries computed from a :class:`~modelcity.settlements.Settlements`.

Every function here returns a plain :class:`pandas.DataFrame` so the numbers can
go into a table or a manuscript without going through a plot first.
"""

from __future__ import annotations

from typing import Iterable, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from .settlements import Settlements

__all__ = [
    "kde",
    "settlement_duration",
    "settlement_events",
    "settlement_churn",
    "size_totals",
    "size_summary",
    "rank_size",
    "urban_share",
]


def kde(
    values: Sequence[float],
    bw: float,
    n: int = 512,
    cut: float = 3.0,
    lower: Optional[float] = None,
    upper: Optional[float] = None,
) -> pd.DataFrame:
    """Kernel density estimate with an absolute bandwidth, as R's ``density``.

    ``scipy.stats.gaussian_kde`` interprets ``bw_method`` as a multiple of the
    sample standard deviation, whereas R's ``bw`` argument is the standard
    deviation of the kernel itself.  Passing ``bw / sd`` reconciles the two, so
    ``kde(x, bw=150)`` matches ``density(x, bw = 150)``.

    The evaluation range also follows R: ``cut`` bandwidths beyond the data.
    """
    array = np.asarray(pd.Series(values).dropna(), dtype=float)
    if array.size == 0:
        return pd.DataFrame({"year": [], "density": []})

    if lower is None:
        lower = float(array.min() - cut * bw)
    if upper is None:
        upper = float(array.max() + cut * bw)
    grid = np.linspace(lower, upper, n)

    spread = float(array.std(ddof=1)) if array.size > 1 else 0.0
    if array.size < 2 or spread == 0.0:
        # A single distinct value has no spread to scale against; place the
        # kernel by hand so the caller still gets a usable curve.
        density = np.exp(-0.5 * ((grid - array.mean()) / bw) ** 2)
        density /= bw * np.sqrt(2.0 * np.pi)
        return pd.DataFrame({"year": grid, "density": density})

    estimator = gaussian_kde(array, bw_method=bw / spread)
    return pd.DataFrame({"year": grid, "density": estimator(grid)})


def settlement_duration(settlements: Settlements) -> pd.DataFrame:
    """Lifespan of each settlement, from first appearance to last.

    For interval data this is the span from the earliest period start to the
    latest period end.  For snapshot data it is the span between the outer
    bounds of the first and last observation.
    """
    occupied = settlements.occupied()
    grouped = occupied.groupby("id", sort=False)
    result = grouped.agg(
        name=("name", "first"),
        group=("group", "first"),
        first_year=("t_start", "min"),
        last_year=("t_end", "max"),
        max_size=("size", "max"),
        n_periods=("size", "size"),
    ).reset_index()
    result["duration"] = result["last_year"] - result["first_year"]
    return result.sort_values("id").reset_index(drop=True)


def settlement_events(settlements: Settlements, bw: float = 150.0, n: int = 512) -> pd.DataFrame:
    """Densities of foundation and abandonment years, on a shared grid."""
    spans = settlement_duration(settlements)
    foundations = spans["first_year"].to_numpy(dtype=float)
    abandonments = spans["last_year"].to_numpy(dtype=float)

    if foundations.size == 0:
        return pd.DataFrame({"year": [], "density": [], "event": []})

    lower = float(min(foundations.min(), abandonments.min()) - 3.0 * bw)
    upper = float(max(foundations.max(), abandonments.max()) + 3.0 * bw)

    found = kde(foundations, bw=bw, n=n, lower=lower, upper=upper)
    found["event"] = "foundations"
    aban = kde(abandonments, bw=bw, n=n, lower=lower, upper=upper)
    aban["event"] = "abandonments"
    return pd.concat([found, aban], ignore_index=True)


def settlement_churn(settlements: Settlements, bw: float = 150.0, n: int = 512) -> pd.DataFrame:
    """Total turnover and its direction over time.

    ``churn`` is foundations plus abandonments, a measure of how much change is
    happening at all.  ``net`` is foundations minus abandonments, positive while
    the settlement system is expanding and negative while it is contracting.
    """
    events = settlement_events(settlements, bw=bw, n=n)
    if events.empty:
        return pd.DataFrame({"year": [], "foundations": [], "abandonments": [], "churn": [], "net": []})

    wide = events.pivot(index="year", columns="event", values="density").reset_index()
    wide.columns.name = None
    wide["churn"] = wide["foundations"] + wide["abandonments"]
    wide["net"] = wide["foundations"] - wide["abandonments"]
    return wide


def _period_frame(settlements: Settlements) -> pd.DataFrame:
    periods = settlements.periods()
    return periods[["t_start", "t_end", "midpoint", "label"]]


def size_totals(settlements: Settlements) -> pd.DataFrame:
    """Total size and settlement count per period."""
    occupied = settlements.occupied()
    totals = (
        occupied.groupby(["t_start", "t_end"], sort=True)
        .agg(total_size=("size", "sum"), n_settlements=("id", "nunique"))
        .reset_index()
    )
    return _period_frame(settlements).merge(totals, on=["t_start", "t_end"], how="left").fillna(
        {"total_size": 0.0, "n_settlements": 0}
    )


def size_summary(settlements: Settlements) -> pd.DataFrame:
    """Median, mean, extremes and the 10-90 percent band per period."""
    occupied = settlements.occupied()
    summary = (
        occupied.groupby(["t_start", "t_end"], sort=True)["size"]
        .agg(
            median="median",
            mean="mean",
            max="max",
            min="min",
            n="size",
            q10=lambda values: values.quantile(0.10),
            q90=lambda values: values.quantile(0.90),
        )
        .reset_index()
    )
    return _period_frame(settlements).merge(summary, on=["t_start", "t_end"], how="left")


def rank_size(settlements: Settlements) -> pd.DataFrame:
    """Rank-size distribution per period, with the Zipf expectation.

    Zipf's law expects the r-th settlement to be 1/r the size of the largest,
    so ``ideal`` is the straight line that a perfectly rank-size ordered system
    would follow on log-log axes.
    """
    occupied = settlements.occupied().copy()
    occupied = occupied.sort_values(["t_start", "t_end", "size"], ascending=[True, True, False])
    occupied["rank"] = occupied.groupby(["t_start", "t_end"]).cumcount() + 1
    largest = occupied.groupby(["t_start", "t_end"])["size"].transform("max")
    occupied["ideal"] = largest / occupied["rank"]

    labels = _period_frame(settlements)
    return occupied.merge(labels, on=["t_start", "t_end"], how="left").reset_index(drop=True)


def urban_share(
    settlements: Settlements,
    thresholds: Optional[Iterable[float]] = None,
) -> pd.DataFrame:
    """Share of the settlement system made up of large centres.

    ``thresholds`` defaults to the pair implied by ``size_type``: 20 and 40
    hectares for areas, 1,000 and 2,000 inhabitants for populations.
    """
    if thresholds is None:
        thresholds = settlements.size.default_thresholds
    thresholds = [float(value) for value in thresholds]

    occupied = settlements.occupied()
    rows = []
    for (start, end), chunk in occupied.groupby(["t_start", "t_end"], sort=True):
        total = float(chunk["size"].sum())
        row = {
            "t_start": start,
            "t_end": end,
            "total_size": total,
            "n_settlements": int(chunk["id"].nunique()),
        }
        for threshold in thresholds:
            large = chunk.loc[chunk["size"] >= threshold, "size"]
            key = "{:g}".format(threshold)
            row["urban_{}".format(key)] = float(large.sum())
            row["n_urban_{}".format(key)] = int(large.count())
            row["percent_{}".format(key)] = (
                float(large.sum()) / total * 100.0 if total > 0 else 0.0
            )
        rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    result = _period_frame(settlements).merge(frame, on=["t_start", "t_end"], how="left")
    result.attrs["thresholds"] = thresholds
    return result
