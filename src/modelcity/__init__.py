"""Unified settlement analysis for long-term urban evolution datasets.

The package standardises any settlement table into one canonical shape, then
computes the same set of summaries and figures from it, whether size is measured
as an area in hectares or as a population count.

    >>> from modelcity import read_yautepec, plot_duration
    >>> sites = read_yautepec("data/examples/yautepec.csv")
    >>> fig = plot_duration(sites)

Datasets other than the two shipped readers go through ``as_settlements``, which
maps your column names onto the canonical roles:

    >>> cities = as_settlements(
    ...     frame, id="Nodes_ID", size="Inhabitants", start="Year",
    ...     size_type="population", size_scale=1000,
    ... )

Years are read in whichever scale the source uses, set by ``time_scale``, and
are written out as years BP.  Use ``set_year_display("BCAD")`` for BC/AD labels
throughout, or ``year_display`` to switch for one figure:

    >>> with year_display("BCAD"):
    ...     fig = plot_events(cities)
"""

from __future__ import annotations

from .metrics import (
    kde,
    rank_size,
    settlement_churn,
    settlement_duration,
    settlement_events,
    size_summary,
    size_totals,
    urban_share,
)
from .plots import (
    plot_churn,
    plot_duration,
    plot_events,
    plot_rank_size,
    plot_size_summary,
    plot_size_totals,
    plot_urban_share,
)
from .readers import read_viabundus, read_yautepec
from .schema import ALL_ROLES, OPTIONAL_ROLES, REQUIRED_ROLES, ColumnRoles, SchemaError
from .settlements import SIZE_TYPES, Settlements, SizeSpec, TimeSpec, as_settlements
from .theme import save_figure
from .timescale import (
    YEAR_DISPLAYS,
    bcad_to_bp,
    bp_to_ce,
    ce_to_bp,
    format_period,
    format_year,
    get_year_display,
    set_year_display,
    to_astronomical,
    to_bp,
    year_display,
)

__version__ = "0.1.0"

__all__ = [
    # construction and validation
    "as_settlements",
    "Settlements",
    "SizeSpec",
    "TimeSpec",
    "SchemaError",
    "ColumnRoles",
    "SIZE_TYPES",
    "REQUIRED_ROLES",
    "OPTIONAL_ROLES",
    "ALL_ROLES",
    # readers
    "read_yautepec",
    "read_viabundus",
    # metrics
    "settlement_duration",
    "settlement_events",
    "settlement_churn",
    "size_totals",
    "size_summary",
    "rank_size",
    "urban_share",
    "kde",
    # plots
    "plot_duration",
    "plot_events",
    "plot_churn",
    "plot_size_totals",
    "plot_size_summary",
    "plot_rank_size",
    "plot_urban_share",
    "save_figure",
    # time helpers
    "bcad_to_bp",
    "bp_to_ce",
    "ce_to_bp",
    "to_astronomical",
    "to_bp",
    "format_year",
    "format_period",
    # year display, BP by default
    "YEAR_DISPLAYS",
    "get_year_display",
    "set_year_display",
    "year_display",
]
