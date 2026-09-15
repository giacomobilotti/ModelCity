"""Shared visual style, so figures from different functions sit together."""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.ticker import FuncFormatter, MaxNLocator

from .timescale import BP_REFERENCE, format_year, get_year_display

__all__ = ["PALETTE", "new_figure", "style_axes", "year_axis", "save_figure"]


PALETTE = {
    "primary": "#8B2323",      # brown4, used for the main density curves
    "secondary": "#CD6839",    # sienna3, used for means and comparisons
    "band": "#CDAA7D",         # burlywood3, used for quantile ribbons
    "foundations": "#7FB8D8",
    "abandonments": "#FF4500",
    "growth": "#006400",
    "decline": "#B22222",
    "accent": "#4682B4",       # steelblue
    "neutral": "#8C8C8C",
    "ideal": "#666666",
    "urban_low": "#3182BD",
    "urban_high": "#DE2D26",
}


def new_figure(
    figsize: Tuple[float, float] = (10.0, 7.5),
    nrows: int = 1,
    ncols: int = 1,
    **kwargs,
):
    """A figure and axes with the package style already applied."""
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize, **kwargs)
    if nrows == 1 and ncols == 1:
        style_axes(axes)
    else:
        for ax in getattr(axes, "flat", [axes]):
            style_axes(ax)
    return fig, axes


def style_axes(ax: Axes) -> Axes:
    """Minimal frame with light horizontal guides, echoing theme_minimal."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#BFBFBF")
    ax.grid(True, axis="both", color="#E6E6E6", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(colors="#4D4D4D", labelsize=9)
    return ax


class _BPLocator(MaxNLocator):
    """Choose ticks that are round numbers of years BP.

    The axis carries astronomical years, so the default locator would pick
    round values on that scale and leave the BP labels reading 3450, 2950,
    2450.  Picking in BP and mapping back keeps the labels legible; the
    mapping is a plain subtraction, so no tick moves off its true position.
    """

    def tick_values(self, vmin, vmax):
        ticks = super().tick_values(BP_REFERENCE - vmax, BP_REFERENCE - vmin)
        return (BP_REFERENCE - ticks)[::-1]


def year_axis(ax: Axes, label: Optional[str] = None) -> Axes:
    """Label the x axis with calendar years in the active display scale.

    Values are astronomical years internally, so the formatter converts on the
    way out; without it a BP axis would read as AD, and a tick at 0 would name
    a year that never existed.  Time still runs left to right under either
    scale, which is why a BP axis counts downwards.
    """
    in_bp = get_year_display() == "BP"
    if in_bp:
        ax.xaxis.set_major_locator(_BPLocator())
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: format_year(value, suffix=False)))
    ax.set_xlabel(label if label is not None else ("Year BP" if in_bp else "Year (BC/AD)"))
    return ax


def save_figure(fig: Figure, path, dpi: int = 150) -> str:
    """Write ``fig`` to ``path``, creating parent directories as needed."""
    from pathlib import Path

    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(destination, dpi=dpi, bbox_inches="tight", facecolor="white")
    return str(destination)
