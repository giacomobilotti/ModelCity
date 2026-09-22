"""The canonical settlement table and its constructor."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .schema import ColumnRoles, SchemaError, resolve_roles
from .timescale import format_period, format_year, snapshot_bounds, to_astronomical, to_bp

__all__ = [
    "Settlements",
    "SizeSpec",
    "TimeSpec",
    "as_settlements",
    "SIZE_TYPES",
]


SIZE_TYPES: Tuple[str, ...] = ("area", "population")

_DEFAULT_UNITS: Dict[str, str] = {"area": "ha", "population": "inhabitants"}

# Smith 2021 treats 20 ha and 40 ha as roughly 1,000 and 2,000 inhabitants, so
# the two sets of defaults describe the same two tiers of settlement.
_DEFAULT_THRESHOLDS: Dict[str, Tuple[float, float]] = {
    "area": (20.0, 40.0),
    "population": (1000.0, 2000.0),
}

_CANONICAL_COLUMNS = ("id", "name", "group", "t_start", "t_end", "size")


@dataclass(frozen=True)
class SizeSpec:
    """What the ``size`` column measures."""

    size_type: str
    unit: str
    scale: float = 1.0

    @property
    def is_area(self) -> bool:
        return self.size_type == "area"

    @property
    def label(self) -> str:
        noun = "area" if self.is_area else "population"
        return "{} [{}]".format(noun.capitalize(), self.unit)

    @property
    def default_thresholds(self) -> Tuple[float, float]:
        return _DEFAULT_THRESHOLDS[self.size_type]


@dataclass(frozen=True)
class TimeSpec:
    """How time was expressed on input and how each row should be read."""

    time_scale: str
    record_type: str

    @property
    def is_snapshot(self) -> bool:
        return self.record_type == "snapshot"


class Settlements:
    """A dataset of settlements observed over time, in canonical form.

    ``data`` always carries the columns ``id``, ``name``, ``group``, ``t_start``,
    ``t_end`` and ``size``.  Years are astronomical calendar years (negative for
    BC) regardless of whether the source counted in BP, and ``size`` is in real
    units after any ``size_scale`` multiplier has been applied.

    The frame is wrapped rather than subclassed: pandas drops unknown attributes
    on most operations, which would silently lose the size and time metadata
    that every downstream function depends on.
    """

    def __init__(self, data: pd.DataFrame, size: SizeSpec, time: TimeSpec, roles: ColumnRoles):
        self.data = data
        self.size = size
        self.time = time
        self.roles = roles

    # -- basic description -------------------------------------------------

    def __len__(self) -> int:
        return len(self.data)

    def __repr__(self) -> str:
        spans = self.data[["t_start", "t_end"]]
        earliest = format_year(spans["t_start"].min())
        latest = format_year(spans["t_end"].max())
        return (
            "<Settlements {n_rows} rows | {n_ids} settlements | {n_periods} periods\n"
            "  size: {size_type} in {unit}\n"
            "  time: {time_scale} input, {record_type} records, {earliest} to {latest}>"
        ).format(
            n_rows=len(self.data),
            n_ids=self.n_settlements,
            n_periods=self.n_periods,
            size_type=self.size.size_type,
            unit=self.size.unit,
            time_scale=self.time.time_scale,
            record_type=self.time.record_type,
            earliest=earliest,
            latest=latest,
        )

    @property
    def n_settlements(self) -> int:
        return int(self.data["id"].nunique())

    @property
    def n_periods(self) -> int:
        return len(self.periods())

    @property
    def is_snapshot(self) -> bool:
        return self.time.is_snapshot

    @property
    def size_label(self) -> str:
        return self.size.label

    @property
    def has_group(self) -> bool:
        return self.data["group"].notna().any()

    def periods(self) -> pd.DataFrame:
        """Distinct periods in chronological order, with display labels."""
        periods = (
            self.data[["t_start", "t_end"]]
            .drop_duplicates()
            .sort_values(["t_start", "t_end"])
            .reset_index(drop=True)
        )
        periods["midpoint"] = (periods["t_start"] + periods["t_end"]) / 2.0
        periods["label"] = [
            format_period(row.t_start, row.t_end) for row in periods.itertuples()
        ]
        return periods

    def occupied(self) -> pd.DataFrame:
        """Rows where the settlement actually existed (``size`` above zero)."""
        return self.data[self.data["size"] > 0].copy()

    def years_bp(self) -> pd.DataFrame:
        """A copy of ``data`` with the year columns restated in years BP.

        BP counts backwards, so ``start_bp`` is the larger number, matching how
        a BP source file is written.  This is for reading and export only: the
        astronomical columns remain the ones every metric is computed from,
        because differences between them are true elapsed spans.
        """
        frame = self.data.copy()
        frame.insert(frame.columns.get_loc("t_start"), "start_bp", to_bp(frame["t_start"]))
        frame.insert(frame.columns.get_loc("t_end"), "end_bp", to_bp(frame["t_end"]))
        if "t_observed" in frame.columns:
            frame["observed_bp"] = to_bp(frame["t_observed"])
        return frame

    def with_data(self, data: pd.DataFrame) -> "Settlements":
        """A copy carrying the same metadata but a different frame."""
        return Settlements(data, self.size, self.time, self.roles)


def _validate_size_type(size_type: Optional[str]) -> str:
    if size_type is None:
        raise SchemaError(
            "size_type is required: pass 'area' for a spatial extent such as "
            "hectares, or 'population' for a headcount.\n"
            "  fix: as_settlements(df, ..., size_type=\"area\")"
        )
    if size_type not in SIZE_TYPES:
        raise SchemaError(
            "size_type must be one of {}, got {!r}.".format(
                " or ".join(repr(value) for value in SIZE_TYPES), size_type
            )
        )
    return size_type


def _require_numeric(frame: pd.DataFrame, column: str, role: str) -> pd.Series:
    values = pd.to_numeric(frame[column], errors="coerce")
    if values.isna().all() and len(values) > 0:
        raise SchemaError(
            "column {!r} (role {!r}) could not be read as numbers.\n"
            "  first values: {}".format(
                column, role, ", ".join(map(str, frame[column].head(3).tolist()))
            )
        )
    return values


def as_settlements(
    data: pd.DataFrame,
    *,
    size_type: Optional[str] = None,
    id: Optional[str] = None,
    size: Optional[str] = None,
    start: Optional[str] = None,
    end: Optional[str] = None,
    name: Optional[str] = None,
    group: Optional[str] = None,
    time_scale: str = "CE",
    size_scale: float = 1.0,
    size_unit: Optional[str] = None,
    drop_missing: bool = True,
) -> Settlements:
    """Standardise ``data`` into a :class:`Settlements` object.

    Parameters
    ----------
    data:
        Any table with one row per settlement per period or observation.
    size_type:
        ``"area"`` or ``"population"``.  Selects axis labels, table units and
        the default urban thresholds.
    id, size, start, end, name, group:
        Column names for each role.  A role is filled automatically when a
        column already carries the canonical name; otherwise it must be given
        here.  ``end`` is omitted for snapshot data.
    time_scale:
        ``"BP"`` if ``start``/``end`` count backwards from 1950, ``"CE"`` if
        they are calendar years.
    size_scale:
        Multiplier bringing ``size`` into real units.  Viabundus stores
        thousands of inhabitants, so it passes ``size_scale=1000``.
    size_unit:
        Overrides the default unit label (``ha`` or ``inhabitants``).
    drop_missing:
        Drop rows with no identifier, size or start year.

    Raises
    ------
    SchemaError
        If a required role cannot be mapped, or if the values in a mapped
        column are not usable.
    """
    if not isinstance(data, pd.DataFrame):
        raise SchemaError(
            "data must be a pandas DataFrame, got {}.".format(type(data).__name__)
        )

    size_type = _validate_size_type(size_type)

    if time_scale not in ("BP", "CE"):
        raise SchemaError(
            "time_scale must be 'BP' or 'CE', got {!r}.".format(time_scale)
        )

    roles = resolve_roles(
        data.columns,
        id=id,
        size=size,
        start=start,
        end=end,
        name=name,
        group=group,
        size_type=size_type,
    )

    frame = pd.DataFrame(index=data.index)
    frame["id"] = data[roles.id].astype(str)
    frame["name"] = (
        data[roles.name].astype(str) if roles.name else frame["id"]
    )
    frame["group"] = data[roles.group] if roles.group else pd.Series(np.nan, index=data.index)

    raw_size = _require_numeric(data, roles.size, "size")
    if (raw_size.dropna() < 0).any():
        offenders = data.loc[raw_size < 0, roles.size].head(3).tolist()
        raise SchemaError(
            "column {!r} (role 'size') contains negative values: {}.\n"
            "  sizes are areas or headcounts and cannot be below zero.".format(
                roles.size, ", ".join(map(str, offenders))
            )
        )
    frame["size"] = raw_size * float(size_scale)

    raw_start = _require_numeric(data, roles.start, "start")
    raw_end = _require_numeric(data, roles.end, "end") if roles.has_end else None

    if time_scale == "BP":
        start_years = to_astronomical(raw_start)
        end_years = to_astronomical(raw_end) if raw_end is not None else None
    else:
        start_years = raw_start.to_numpy(dtype=float)
        end_years = raw_end.to_numpy(dtype=float) if raw_end is not None else None

    if end_years is None:
        record_type = "snapshot"
        bounds = snapshot_bounds(start_years[np.isfinite(start_years)])
        lower = np.array([bounds.get(year, (year, year))[0] for year in start_years])
        upper = np.array([bounds.get(year, (year, year))[1] for year in start_years])
        frame["t_start"] = lower
        frame["t_end"] = upper
        frame["t_observed"] = start_years
    else:
        record_type = "interval"
        frame["t_start"] = start_years
        frame["t_end"] = end_years
        frame["t_observed"] = (start_years + end_years) / 2.0

    # Preserve anything else the caller may want later, such as coordinates.
    used = {roles.id, roles.size, roles.start}
    for optional in (roles.end, roles.name, roles.group):
        if optional:
            used.add(optional)
    for column in data.columns:
        if column not in used and column not in frame.columns:
            frame[column] = data[column]

    if drop_missing:
        frame = frame.dropna(subset=["id", "size", "t_start", "t_end"])

    reversed_rows = frame["t_start"] > frame["t_end"]
    if reversed_rows.any():
        example = frame.loc[reversed_rows, ["id", "t_start", "t_end"]].head(3)
        raise SchemaError(
            "{} row(s) end before they start once converted to calendar years.\n"
            "  with time_scale='BP' the start year is the larger number, "
            "because BP counts backwards.\n"
            "  first offenders:\n{}".format(int(reversed_rows.sum()), example.to_string(index=False))
        )

    ordered = [column for column in _CANONICAL_COLUMNS if column in frame.columns]
    remainder = [column for column in frame.columns if column not in ordered]
    frame = frame[ordered + remainder].reset_index(drop=True)

    spec = SizeSpec(
        size_type=size_type,
        unit=size_unit or _DEFAULT_UNITS[size_type],
        scale=float(size_scale),
    )
    time = TimeSpec(time_scale=time_scale, record_type=record_type)
    return Settlements(frame, spec, time, roles)
