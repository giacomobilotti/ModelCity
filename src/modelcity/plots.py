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

from . import metrics
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
