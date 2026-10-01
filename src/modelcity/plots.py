"""The core figures, each returning a :class:`matplotlib.figure.Figure`.

Axis labels and default thresholds come from the object's metadata, so the same
call produces "Total occupied area [ha]" for an area dataset and "Total
population [inhabitants]" for a population one.
"""

from __future__ import annotations

import math
from typing import Iterable, Optional, Tuple

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from . import metrics, sampling
from .settlements import Settlements
from .theme import PALETTE, new_figure, style_axes, year_axis

__all__ = [
    "plot_duration",
    "plot_events",
    "plot_churn",
    "plot_size_totals",
    "plot_size_summary",
    "plot_rank_size",
    "plot_urban_share",
    "plot_site_growth",
    "plot_trajectories",
    "plot_period_growth",
    "plot_persistence",
]


def _noun(settlements: Settlements) -> str:
    return "site" if settlements.size.is_area else "city"


def _total_label(settlements: Settlements) -> str:
    if settlements.size.is_area:
        return "Total occupied area [{}]".format(settlements.size.unit)
    return "Total population [{}]".format(settlements.size.unit)


def _find_peaks(x: np.ndarray, y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Local maxima, where the slope turns from rising to falling."""
    if y.size < 3:
        return np.array([]), np.array([])
    turning = np.diff(np.sign(np.diff(y)))
    indices = np.where(turning == -2)[0] + 1
    return x[indices], y[indices]


def plot_duration(
    settlements: Settlements,
    bw: float = 150.0,
    mark_peaks: bool = True,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """Distribution of how long settlements lasted."""
    spans = metrics.settlement_duration(settlements)
    density = metrics.kde(spans["duration"], bw=bw)

    fig, ax = new_figure(figsize=figsize)
    ax.plot(density["year"], density["density"], color=PALETTE["primary"], linewidth=2)

    if mark_peaks and not density.empty:
        peak_x, peak_y = _find_peaks(
            density["year"].to_numpy(), density["density"].to_numpy()
        )
        if peak_x.size:
            ax.plot(peak_x, peak_y, "o", color="#8B0000", markersize=5)
            for x_value, y_value in zip(peak_x, peak_y):
                ax.annotate(
                    "{:.0f}".format(x_value),
                    (x_value, y_value),
                    textcoords="offset points",
                    xytext=(0, 7),
                    ha="center",
                    fontsize=8,
                    color="#8B0000",
                )

    median = float(spans["duration"].median()) if len(spans) else float("nan")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Duration (years)")
    ax.set_ylabel("Density")
    ax.set_title(
        "{} duration (median {:.0f} years, n = {})".format(
            _noun(settlements).capitalize(), median, len(spans)
        )
    )
    fig.tight_layout()
    return fig


def plot_events(
    settlements: Settlements,
    bw: float = 150.0,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """When settlements were founded and when they were abandoned."""
    events = metrics.settlement_events(settlements, bw=bw)

    fig, ax = new_figure(figsize=figsize)
    for event, colour in (
        ("foundations", PALETTE["foundations"]),
        ("abandonments", PALETTE["abandonments"]),
    ):
        subset = events[events["event"] == event]
        ax.plot(
            subset["year"],
            subset["density"],
            color=colour,
            linewidth=1.8,
            label="{} {}".format(_noun(settlements).capitalize(), event),
        )
        ax.fill_between(subset["year"], subset["density"], color=colour, alpha=0.15)

    ax.set_ylim(bottom=0)
    ax.set_ylabel("Density")
    year_axis(ax)
    ax.legend(frameon=False)
    ax.set_title("Foundations and abandonments (bandwidth {:g} years)".format(bw))
    fig.tight_layout()
    return fig


def plot_churn(
    settlements: Settlements,
    bw: float = 150.0,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """Total turnover, and whether foundations or abandonments dominate."""
    churn = metrics.settlement_churn(settlements, bw=bw)

    fig, ax = new_figure(figsize=figsize)
    years = churn["year"].to_numpy()
    net = churn["net"].to_numpy()

    ax.fill_between(
        years, 0, np.where(net > 0, net, 0), color=PALETTE["growth"], alpha=0.25,
        label="Net growth",
    )
    ax.fill_between(
        years, np.where(net < 0, net, 0), 0, color=PALETTE["decline"], alpha=0.25,
        label="Net decline",
    )
    ax.axhline(0, linestyle="--", color="#808080", linewidth=1)
    ax.plot(years, churn["churn"], color=PALETTE["accent"], linewidth=1.6, label="Churn (total change)")
    ax.plot(years, net, color=PALETTE["neutral"], linewidth=1.4, label="Net (foundations - abandonments)")

    ax.set_ylabel("Density")
    year_axis(ax)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Settlement turnover")
    fig.tight_layout()
    return fig


def plot_size_totals(
    settlements: Settlements,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """Total size per period, with the settlement count on a second axis."""
    totals = metrics.size_totals(settlements)

    fig, ax = new_figure(figsize=figsize)
    for row in totals.itertuples():
        ax.hlines(
            row.total_size,
            row.t_start,
            row.t_end,
            color="#404040",
            linewidth=6,
            alpha=0.85,
        )

    counts = ax.twinx()
    counts.plot(
        totals["midpoint"],
        totals["n_settlements"],
        marker="^",
        markersize=8,
        markerfacecolor="#B0D4E3",
        markeredgecolor="black",
        linestyle="none",
    )
    counts.set_ylabel("Number of {}s".format(_noun(settlements)))
    counts.set_ylim(bottom=0)
    counts.grid(False)
    for side in ("top", "left"):
        counts.spines[side].set_visible(False)

    ax.set_ylim(bottom=0)
    ax.set_ylabel(_total_label(settlements))
    year_axis(ax)
    ax.set_title("{} totals and counts over time".format(_noun(settlements).capitalize()))
    fig.tight_layout()
    return fig


def plot_size_summary(
    settlements: Settlements,
    figsize: Tuple[float, float] = (10.0, 7.5),
    log_scale: bool = False,
) -> Figure:
    """Median and mean size per period against the 10-90 percent band."""
    summary = metrics.size_summary(settlements)

    fig, ax = new_figure(figsize=figsize)
    for row in summary.itertuples():
        ax.fill_between(
            [row.t_start, row.t_end],
            [row.q10, row.q10],
            [row.q90, row.q90],
            color=PALETTE["band"],
            alpha=0.35,
            linewidth=0,
        )
        ax.hlines(row.median, row.t_start, row.t_end, color=PALETTE["primary"], linewidth=2.2)
        ax.hlines(
            row.mean, row.t_start, row.t_end,
            color=PALETTE["secondary"], linewidth=1.6, linestyle="--",
        )

    ax.plot(
        summary["midpoint"], summary["max"],
        marker="o", markersize=5, markerfacecolor="white",
        markeredgecolor="black", color="#4D4D4D", linewidth=1,
        label="Largest {}".format(_noun(settlements)),
    )

    if log_scale:
        ax.set_yscale("log")

    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D([], [], color=PALETTE["primary"], linewidth=2.2, label="Median"),
            Line2D([], [], color=PALETTE["secondary"], linewidth=1.6, linestyle="--", label="Mean"),
            Line2D([], [], color=PALETTE["band"], linewidth=8, alpha=0.5, label="10-90 percent"),
            Line2D(
                [], [], color="#4D4D4D", marker="o", markerfacecolor="white",
                markeredgecolor="black", label="Largest {}".format(_noun(settlements)),
            ),
        ],
        frameon=False,
        loc="upper left",
    )
    ax.set_ylabel(settlements.size.label)
    year_axis(ax)
    ax.set_title("{} size distribution over time".format(_noun(settlements).capitalize()))
    fig.tight_layout()
    return fig


def plot_rank_size(
    settlements: Settlements,
    ncol: int = 4,
    figsize: Optional[Tuple[float, float]] = None,
) -> Figure:
    """Rank-size curves per period against the Zipf expectation."""
    ranked = metrics.rank_size(settlements)
    periods = settlements.periods()
    n_periods = len(periods)
    if n_periods == 0:
        raise ValueError("no periods to plot")

    ncol = max(1, min(ncol, n_periods))
    nrow = math.ceil(n_periods / ncol)
    if figsize is None:
        figsize = (3.2 * ncol, 2.8 * nrow)

    fig, axes = new_figure(figsize=figsize, nrows=nrow, ncols=ncol, squeeze=False)
    flat = axes.flat

    for index, period in enumerate(periods.itertuples()):
        ax = flat[index]
        subset = ranked[
            (ranked["t_start"] == period.t_start) & (ranked["t_end"] == period.t_end)
        ]
        # Markers matter here: a period with a single settlement would draw an
        # empty panel if it were rendered as a line alone.
        ax.plot(
            subset["rank"], subset["size"],
            color=PALETTE["accent"], linewidth=1.4, marker="o", markersize=2.5,
        )
        ax.plot(subset["rank"], subset["ideal"], color=PALETTE["ideal"], linewidth=1, linestyle="--")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(period.label, fontsize=9)
        ax.tick_params(labelsize=8)

    for leftover in range(n_periods, nrow * ncol):
        flat[leftover].set_visible(False)

    fig.supxlabel("Rank (log scale)")
    fig.supylabel("{} (log scale)".format(settlements.size.label))
    fig.suptitle("Rank-size distribution by period")
    fig.tight_layout()
    return fig


def plot_urban_share(
    settlements: Settlements,
    thresholds: Optional[Iterable[float]] = None,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """Absolute and relative weight of large centres over time."""
    shares = metrics.urban_share(settlements, thresholds=thresholds)
    used = shares.attrs.get("thresholds", settlements.size.default_thresholds)
    colours = [PALETTE["urban_low"], PALETTE["urban_high"], PALETTE["growth"]]

    fig, ax = new_figure(figsize=figsize)
    percent = ax.twinx()

    for index, threshold in enumerate(used):
        key = "{:g}".format(threshold)
        colour = colours[index % len(colours)]
        label = ">= {:g} {}".format(threshold, settlements.size.unit)
        ax.plot(
            shares["midpoint"], shares["urban_{}".format(key)],
            color=colour, linewidth=1.8, marker="o", markersize=4, label=label,
        )
        percent.plot(
            shares["midpoint"], shares["percent_{}".format(key)],
            color=colour, linewidth=1.4, linestyle="--", alpha=0.8,
        )

    ax.set_ylim(bottom=0)
    percent.set_ylim(0, 100)
    percent.set_ylabel("Percentage of total (dashed)")
    percent.grid(False)
    for side in ("top", "left"):
        percent.spines[side].set_visible(False)

    ax.set_ylabel(_total_label(settlements))
    year_axis(ax)
    ax.legend(frameon=False, loc="upper left", title="Threshold")
    ax.set_title("Weight of large centres over time")
    fig.tight_layout()
    return fig


# -- chronological sampling --------------------------------------------------
#
# These take the frame returned by ``sampling.sample_trajectories`` rather than
# a Settlements object.  Its ``attrs`` carry the size labels, which pandas
# drops on some operations, so every label has a neutral fallback.

_SKYBLUE = "#87CEEB"

# Trajectory column and the site_summary column drawn as its median.
_TRAJECTORY_VALUES = {"growth_cagr": "median", "growth_linear": "med_lin", "size": "size"}


def _trajectory_noun(trajectories: pd.DataFrame) -> str:
    return "city" if trajectories.attrs.get("size_type") == "population" else "site"


def _one_site(trajectories: pd.DataFrame, site_id) -> pd.DataFrame:
    site = trajectories[trajectories["id"].astype(str) == str(site_id)]
    if site.empty:
        raise ValueError("no trajectories for id {!r}.".format(site_id))
    return site


def _site_span(summary: pd.DataFrame) -> Tuple[float, float]:
    return float(summary["t_start"].min()), float(summary["t_end"].max())


def _r_limits(low: float, high: float) -> Tuple[float, float]:
    """Widen limits by 4 percent each side, as R's default axis style does."""
    pad = 0.04 * (high - low)
    return low - pad, high + pad


def plot_site_growth(
    trajectories: pd.DataFrame,
    site_id,
    legend: bool = True,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """R's ``plot_site_cagr()``: median and mean CAGR of one settlement, with 50 and 95 percent bands.

    Points sit at the median sampled date of each phase from
    :func:`~modelcity.sampling.site_summary`.  As in R, the bands and the y
    range leave out the earliest phase, whose median is missing because the
    foundation falls in it, so a single-phase settlement cannot be drawn.
    """
    summary = sampling.site_summary(_one_site(trajectories, site_id))
    if len(summary) < 2:
        raise ValueError("id {!r} has a single phase, which leaves no growth band to plot.".format(site_id))
    # site_summary lists the latest phase first, as R does, so the earliest is last.
    body = summary.iloc[:-1].sort_values("year_median")
    ordered = summary.sort_values("year_median")

    fig, ax = new_figure(figsize=figsize)
    ax.fill_between(
        body["year_median"], body["q025"], body["q975"],
        color=PALETTE["accent"], alpha=0.2, linewidth=0, label="95% CI",
    )
    ax.fill_between(
        body["year_median"], body["q25"], body["q75"],
        color=PALETTE["accent"], alpha=0.3, linewidth=0, label="50% CI",
    )
    ax.plot(ordered["year_median"], ordered["median"], color=PALETTE["accent"], linewidth=2, marker="o", label="Median")
    ax.plot(ordered["year_median"], ordered["mean"], color=_SKYBLUE, linewidth=2, linestyle="--", label="Mean")
    ax.axhline(0, linestyle="--", color="#1A1A1A", linewidth=1)

    low, high = np.nanmin(body["q025"]), np.nanmax(body["q975"])
    if np.isfinite(low) and np.isfinite(high) and high > low:
        ax.set_ylim(*_r_limits(low, high))
    ax.set_xlim(*_r_limits(*_site_span(summary)))
    ax.set_ylabel("CAGR")
    year_axis(ax)
    if legend:
        ax.legend(frameon=False, loc="upper left")
    ax.set_title("Growth Rate for Site {}".format(site_id))
    fig.tight_layout()
    return fig


def plot_trajectories(
    trajectories: pd.DataFrame,
    site_id,
    value: str = "growth_cagr",
    n: int = 50,
    seed=1234,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """R's ``plot_spaghetti()``: ``n`` random draws of one settlement, with the median on top.

    ``value`` is ``"growth_cagr"``, ``"growth_linear"`` or ``"size"`` (R's
    ``"area"`` is accepted too).  The median line is the matching
    :func:`~modelcity.sampling.site_summary` column, and dashed vertical lines
    mark each phase start.  The default y range is R's: the data range
    widened by a tenth of each limit.  ``xlim`` takes astronomical years.
    """
    value = "size" if value == "area" else value
    if value not in _TRAJECTORY_VALUES:
        raise ValueError(
            "value must be 'growth_cagr', 'growth_linear' or 'size', got {!r}.".format(value)
        )
    site = _one_site(trajectories, site_id)
    summary = sampling.site_summary(site)

    draws = np.unique(site["iteration"])
    rng = np.random.default_rng(seed)
    chosen = rng.choice(draws, size=min(n, draws.size), replace=False)

    fig, ax = new_figure(figsize=figsize)
    for _, draw in site[site["iteration"].isin(chosen)].groupby("iteration"):
        draw = draw.sort_values("year", kind="stable")
        ax.plot(draw["year"], draw[value], color="#000000", alpha=0.2, linewidth=0.8)

    if value != "size":
        ax.axhline(0, linestyle="--", color="#808080", linewidth=1)
    ordered = summary.sort_values("year_median")
    ax.plot(ordered["year_median"], ordered[_TRAJECTORY_VALUES[value]], color=_SKYBLUE, linewidth=3)
    for start in summary["t_start"]:
        ax.axvline(start, linestyle="--", color="#000000", linewidth=0.75)

    if ylim is None:
        low, high = np.nanmin(site[value]), np.nanmax(site[value])
        ylim = (low - 0.1 * low, high + 0.1 * high)
    ax.set_xlim(*_r_limits(*(xlim if xlim is not None else _site_span(summary))))
    ax.set_ylim(*_r_limits(*ylim))

    label = trajectories.attrs.get("size_type", "size") if value == "size" else value
    ax.set_ylabel(label)
    year_axis(ax)
    ax.set_title(
        "Site {} {} trajectories ({} samples)".format(site_id, label.split("_")[0], len(chosen))
    )
    fig.tight_layout()
    return fig


def plot_period_growth(
    trajectories: pd.DataFrame,
    figsize: Tuple[float, float] = (10.0, 7.5),
) -> Figure:
    """Median CAGR per phase from :func:`~modelcity.sampling.period_summary`, with its 95 percent band."""
    growth = sampling.period_summary(trajectories).sort_values("year_median")

    fig, ax = new_figure(figsize=figsize)
    ax.fill_between(
        growth["year_median"], growth["q025"], growth["q975"],
        color=PALETTE["accent"], alpha=0.2, linewidth=0, label="95% CI",
    )
    ax.plot(growth["year_median"], growth["median"], color=PALETTE["accent"], linewidth=2, marker="o", label="Median")
    ax.axhline(0, linestyle="--", color="#1A1A1A", linewidth=1)
    ax.set_ylabel("CAGR")
    year_axis(ax)
    ax.legend(frameon=False, loc="upper left")
    ax.set_title("Median CAGR by period")
    fig.tight_layout()
    return fig


def plot_persistence(
    trajectories: pd.DataFrame,
    terminus: float,
    normalised: bool = False,
    time_scale: Optional[str] = None,
    figsize: Tuple[float, float] = (11.0, 7.0),
) -> Figure:
    """Median persistence to ``terminus`` of each settlement, with its 95 percent interval.

    ``terminus`` and ``time_scale`` are as in
    :func:`~modelcity.sampling.persistence_samples`.  One bar per settlement,
    so subset ``trajectories`` first for large datasets.
    """
    summary = sampling.persistence_summary(trajectories, terminus, time_scale)
    suffix = "_n" if normalised else ""
    summary = summary.sort_values("median" + suffix, ascending=False).reset_index(drop=True)

    median = summary["median" + suffix].to_numpy()
    errors = np.vstack([median - summary["q025" + suffix], summary["q975" + suffix] - median])
    positions = np.arange(len(summary))

    fig, ax = new_figure(figsize=figsize)
    ax.bar(positions, median, color=PALETTE["accent"], yerr=errors, capsize=3, ecolor="#404040")
    ax.set_xticks(positions)
    ax.set_xticklabels(summary["id"], rotation=90, fontsize=8)
    ax.set_xlabel(_trajectory_noun(trajectories).capitalize())
    ax.set_ylabel("Share of available time" if normalised else "Median persistence (years)")
    ax.set_ylim(bottom=0)
    ax.set_title("{} persistence to terminus {:g}".format(_trajectory_noun(trajectories).capitalize(), terminus))
    fig.tight_layout()
    return fig
