"""Loaders for the two datasets in this project.

These are thin conveniences: they only supply the role mapping and metadata that
:func:`~modelcity.settlements.as_settlements` would otherwise need spelled out.
Any other dataset can be used directly by calling ``as_settlements`` yourself.
"""

from __future__ import annotations

import csv
import re
import warnings
from io import StringIO
from typing import Optional, Tuple

import numpy as np
import pandas as pd

from .settlements import Settlements, as_settlements

__all__ = [
    "read_yautepec",
    "read_viabundus",
    "read_r_csv",
    "parse_r_point",
    "assign_regions",
]


#: Viabundus stores geometry as the R printout of a coordinate pair.
_POINT_PATTERN = re.compile(r"c\(\s*(-?[0-9.]+)\s*,\s*(-?[0-9.]+)\s*\)")

#: The Viabundus geometry column is in ETRS89 / LAEA Europe.
VIABUNDUS_CRS = 3035


def parse_r_point(values: pd.Series) -> pd.DataFrame:
    """Extract x and y from strings shaped like ``c(4316502.65, 3158017.76)``."""
    extracted = values.astype(str).str.extract(_POINT_PATTERN)
    return pd.DataFrame(
        {
            "x": pd.to_numeric(extracted[0], errors="coerce"),
            "y": pd.to_numeric(extracted[1], errors="coerce"),
        },
        index=values.index,
    )


def _is_split_point(before: str, after: str) -> bool:
    """True when two adjacent fields are the halves of one ``c(x, y)``."""
    return before.strip().startswith("c(") and after.strip().endswith(")")


def read_r_csv(path) -> pd.DataFrame:
    """Read a CSV exported from R whose last column is an unquoted ``c(x, y)``.

    The comma inside that coordinate pair splits it across two fields, leaving
    each row with one more field than the header has names.  pandas resolves the
    mismatch by promoting the first column to the index, which shifts every
    other column one place to the left.  That is easy to miss, because the
    numbers still look plausible; they just belong to the neighbouring column.

    Records are read with the :mod:`csv` module rather than repaired as text,
    because description fields in these exports contain embedded newlines, so a
    physical line is not the same thing as a row.
    """
    with open(path, "r", encoding="utf-8", newline="") as handle:
        rows = list(csv.reader(handle))

    if not rows:
        return pd.DataFrame()

    header = rows[0]
    width = len(header)

    repaired = []
    for row in rows[1:]:
        while len(row) > width and _is_split_point(row[-2], row[-1]):
            row = row[:-2] + [", ".join(field.strip() for field in row[-2:])]
        repaired.append(row)

    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(header)
    writer.writerows(repaired)
    buffer.seek(0)

    frame = pd.read_csv(buffer, low_memory=False)

    # R writes row numbers into an unnamed leading column.
    if len(frame.columns) and str(frame.columns[0]).startswith("Unnamed:"):
        frame = frame.drop(columns=frame.columns[0])
    return frame


def assign_regions(lon: pd.Series, lat: pd.Series) -> pd.Series:
    """Approximate ethno-linguistic regions from position.

    A stand-in for the project's ``linguistic_areas.gpkg``, which is not part of
    this repository.  Boundaries are rough and meant for grouping figures, not
    for any claim about where languages were actually spoken.
    """
    region = pd.Series("germanic", index=lon.index, dtype=object)
    region[(lon >= 14.5) & (lat < 55.5)] = "slavic"
    region[(lon >= 7.5) & (lon < 16.5) & (lat >= 54.0)] = "danish"
    region[(lon < 7.2) & (lat >= 50.0) & (lat < 54.5)] = "dutch"
    region[lon.isna() | lat.isna()] = np.nan
    return region


def _to_lonlat(x: pd.Series, y: pd.Series, crs: int) -> Optional[Tuple[pd.Series, pd.Series]]:
    try:
        from pyproj import Transformer
    except ImportError:
        return None
    transformer = Transformer.from_crs(crs, 4326, always_xy=True)
    lon, lat = transformer.transform(x.to_numpy(), y.to_numpy())
    return pd.Series(lon, index=x.index), pd.Series(lat, index=y.index)


def read_yautepec(path, **kwargs) -> Settlements:
    """Load the Yautepec survey: site areas in hectares, dated in years BP."""
    frame = pd.read_csv(path)
    options = dict(
        id="Sitio",
        size="Area",
        start="start",
        end="end",
        size_type="area",
        time_scale="BP",
    )
    options.update(kwargs)
    return as_settlements(frame, **options)


def read_viabundus(path, assign_region: bool = True, **kwargs) -> Settlements:
    """Load the Viabundus cities: inhabitants in thousands, snapshot years CE.

    Population is stored in thousands, so ``size_scale=1000`` brings it into
    real headcounts and makes the population thresholds directly meaningful.
    """
    frame = read_r_csv(path)

    if "geometry" in frame.columns:
        points = parse_r_point(frame["geometry"])
        frame["x"] = points["x"]
        frame["y"] = points["y"]

        if assign_region:
            lonlat = _to_lonlat(frame["x"], frame["y"], VIABUNDUS_CRS)
            if lonlat is None:
                warnings.warn(
                    "pyproj is not installed, so regions were not assigned. "
                    "Install the geo extra with: pip install 'modelcity[geo]'",
                    stacklevel=2,
                )
            else:
                frame["lon"], frame["lat"] = lonlat
                frame["region"] = assign_regions(frame["lon"], frame["lat"])

    options = dict(
        id="Nodes_ID",
        size="Inhabitants",
        start="Year",
        name="Name",
        size_type="population",
        time_scale="CE",
        size_scale=1000.0,
    )
    if "region" in frame.columns:
        options["group"] = "region"
    options.update(kwargs)
    return as_settlements(frame, **options)
