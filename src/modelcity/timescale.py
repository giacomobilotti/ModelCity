"""Conversion between years BP and calendar years, plus snapshot binning.

Three numberings are in play, and keeping them apart matters:

*BP* counts backwards from 1950 CE.  It is the default for reading and for
display, because both source datasets and the wider literature are quoted this
way.

*Astronomical* numbering has a year zero, so ``1 BC`` is ``0`` and ``2 BC`` is
``-1``.  Differences between astronomical years are true elapsed spans, which is
what durations, midpoints and density grids need, so this is what the package
stores internally.

*Historical* BC/AD numbering has no year zero: ``1 BC`` is followed directly by
``1 AD``.  It is only ever produced for labels.

Astronomical years relate to BP by a single subtraction, ``bp = 1950 - year``,
with no special case at the boundary.  That is the whole reason the internal
store is astronomical: the year-zero correction happens once, when a label is
formatted, instead of inside every calculation.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Dict, Iterator, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

__all__ = [
    "BP_REFERENCE",
    "YEAR_DISPLAYS",
    "bcad_to_bp",
    "bp_to_ce",
    "ce_to_bp",
    "to_astronomical",
    "to_bp",
    "format_year",
    "format_period",
    "get_year_display",
    "set_year_display",
    "year_display",
    "snapshot_bounds",
]

#: Years BP are counted back from 1950 CE by convention.
BP_REFERENCE = 1950

#: Recognised display scales for axis ticks and period labels.
YEAR_DISPLAYS: Tuple[str, ...] = ("BP", "BCAD")

ArrayLike = Union[float, int, Sequence[float], np.ndarray, pd.Series]

_year_display = "BP"


# -- the R conversion, ported ----------------------------------------------

def bcad_to_bp(x: ArrayLike, bc_to_bp: bool = True) -> np.ndarray:
    """Convert between BC/AD and BP, following the project's R implementation.

    This is a direct port of the ``bcad_to_bp()`` used alongside the R scripts,
    which merges ``rcarbon::BCADtoBP`` and ``rcarbon::BPtoBCAD`` into one
    function and drops the upper guard on AD dates.  Both directions keep the
    no-year-zero rule, so ``1 BC`` is 1950 BP and ``1 AD`` is 1949 BP.

    Parameters
    ----------
    x:
        Years to convert.  ``NaN`` passes through untouched, as ``NA`` does in R.
    bc_to_bp:
        ``True`` converts BC/AD to BP, ``False`` converts BP to BC/AD.

    Raises
    ------
    ValueError
        On year 0 when converting from BC/AD, since it does not exist, and on
        negative BP when converting the other way, since post-bomb dates are
        out of scope.

    Notes
    -----
    The R version reaches AD dates after 1950 through ``abs()``, which returns
    a positive BP value where a negative one is meant.  That is reproduced here
    for parity; neither shipped dataset goes near the boundary.
    """
    values = np.atleast_1d(np.asarray(x, dtype=float))
    known = ~np.isnan(values)
    result = np.full(values.shape, np.nan, dtype=float)

    if bc_to_bp:
        if np.any(values[known] == 0):
            raise ValueError("0 BC/AD is not a valid year.")
        ad = known & (values > 0)
        bc = known & (values < 0)
        result[ad] = np.abs(values[ad] - BP_REFERENCE)
        result[bc] = np.abs(values[bc] - (BP_REFERENCE - 1))
    else:
        if np.any(values[known] < 0):
            raise ValueError("Post-bomb dates (<0 BP) are not currently supported.")
        recent = known & (values < BP_REFERENCE)
        ancient = known & (values >= BP_REFERENCE)
        result[recent] = BP_REFERENCE - values[recent]
        result[ancient] = (BP_REFERENCE - 1) - values[ancient]

    return result if np.ndim(x) else result[0]


def ce_to_bp(ce: ArrayLike) -> np.ndarray:
    """Historical BC/AD years to years BP.  Calls :func:`bcad_to_bp`."""
    return bcad_to_bp(ce, bc_to_bp=True)


def bp_to_ce(bp: ArrayLike) -> np.ndarray:
    """Years BP to historical BC/AD years.  Calls :func:`bcad_to_bp`.

    Intended for labelling; prefer :func:`to_astronomical` when the result will
    be used in a calculation.
    """
    return bcad_to_bp(bp, bc_to_bp=False)


# -- the internal linear scale ----------------------------------------------

def to_astronomical(bp: ArrayLike) -> np.ndarray:
    """Years BP to astronomical calendar years, keeping a year zero.

    Use this for arithmetic.  ``to_astronomical(2049) - to_astronomical(1751)``
    gives 298 years, matching the plain difference of the BP values.
    """
    return BP_REFERENCE - np.asarray(bp, dtype=float)


def to_bp(year: ArrayLike) -> np.ndarray:
    """Astronomical calendar years back to years BP, inverting :func:`to_astronomical`.

    Unlike :func:`ce_to_bp` this takes astronomical years, so it needs no
    year-zero correction and stays exact for fractional values such as period
    midpoints.
    """
    return BP_REFERENCE - np.asarray(year, dtype=float)


# -- display scale ----------------------------------------------------------

def get_year_display() -> str:
    """The scale currently used for axis ticks and period labels."""
    return _year_display


def set_year_display(display: str) -> str:
    """Set the display scale to ``"BP"`` or ``"BCAD"``, returning the previous one.

    Only labels are affected.  Stored years stay astronomical either way, so
    switching cannot change any computed number.
    """
    global _year_display
    if display not in YEAR_DISPLAYS:
        raise ValueError(
            "year display must be one of {}, got {!r}.".format(
                " or ".join(repr(value) for value in YEAR_DISPLAYS), display
            )
        )
    previous = _year_display
    _year_display = display
    return previous


@contextmanager
def year_display(display: str) -> Iterator[str]:
    """Use a display scale for one block, then restore the previous one.

        >>> with year_display("BCAD"):
        ...     fig = plot_events(cities)
    """
    previous = set_year_display(display)
    try:
        yield display
    finally:
        set_year_display(previous)


def _resolve_display(display: Optional[str]) -> str:
    return _year_display if display is None else display


def format_year(
    year: float,
    suffix: bool = True,
    display: Optional[str] = None,
) -> str:
    """Format an astronomical year for reading, in the active display scale.

    ``display`` overrides the package default for this call alone.
    """
    if not np.isfinite(year):
        return ""
    value = float(year)

    if _resolve_display(display) == "BP":
        years_bp = BP_REFERENCE - value
        return "{:g} BP".format(years_bp) if suffix else "{:g}".format(years_bp)

    if value < 1:
        # Astronomical 0 is 1 BC, -1 is 2 BC, and so on.
        label = "{:g}".format(abs(value - 1))
        return "{} BC".format(label) if suffix else "-{}".format(label)
    return "{:g} AD".format(value) if suffix else "{:g}".format(value)


def format_period(start: float, end: float, display: Optional[str] = None) -> str:
    """Label a period given its astronomical bounds.

    The era is written once when both ends share it, so a period reads
    ``"3449-3050 BP"`` rather than repeating the suffix.
    """
    scale = _resolve_display(display)
    first = format_year(start, display=scale)
    last = format_year(end, display=scale)
    if not first or not last:
        return first or last

    era = first.rsplit(" ", 1)[-1]
    if era == last.rsplit(" ", 1)[-1]:
        return "{}-{}".format(first.rsplit(" ", 1)[0], last)
    return "{} to {}".format(first, last)


def snapshot_bounds(years: Sequence[float]) -> Dict[float, Tuple[float, float]]:
    """Turn observation years into contiguous intervals via midpoints.

    Snapshot datasets record a settlement at discrete moments (Viabundus has
    1300, 1400, 1500, 1550, 1600, 1650) rather than over explicit periods.  Each
    observation is taken to represent the span halfway to its neighbours; the
    first and last are extended outwards by half of their single neighbouring
    gap so no part of the timeline is unassigned.
    """
    unique = sorted({float(year) for year in years if np.isfinite(year)})
    if not unique:
        return {}
    if len(unique) == 1:
        only = unique[0]
        return {only: (only, only)}

    midpoints = [
        (unique[index] + unique[index + 1]) / 2.0 for index in range(len(unique) - 1)
    ]

    bounds: Dict[float, Tuple[float, float]] = {}
    for index, year in enumerate(unique):
        if index == 0:
            lower = year - (midpoints[0] - year)
        else:
            lower = midpoints[index - 1]
        if index == len(unique) - 1:
            upper = year + (year - midpoints[-1])
        else:
            upper = midpoints[index]
        bounds[year] = (lower, upper)
    return bounds
