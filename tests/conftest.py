from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

from modelcity.timescale import get_year_display, set_year_display

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "examples"
SCRIPTS = ROOT / "scripts"

# Dataset loaders live in scripts/run_unified.py, not in the modelcity package.
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


@pytest.fixture(autouse=True)
def isolate_year_display():
    """The display scale is module-level state, so no test may leak its choice."""
    previous = get_year_display()
    yield
    set_year_display(previous)


@pytest.fixture
def interval_frame() -> pd.DataFrame:
    """A tiny Yautepec-shaped dataset: areas in hectares, dated in years BP."""
    return pd.DataFrame(
        {
            "Sitio": ["A", "A", "B", "C"],
            "Area": [1.0, 3.0, 50.0, 0.5],
            "start": [3449, 3049, 3049, 1750],
            "end": [3050, 2450, 2450, 1651],
        }
    )


@pytest.fixture
def snapshot_frame() -> pd.DataFrame:
    """A tiny Viabundus-shaped dataset: populations recorded at set years."""
    return pd.DataFrame(
        {
            "Nodes_ID": [1, 1, 1, 2, 2, 2],
            "Name": ["Alpha"] * 3 + ["Beta"] * 3,
            "Year": [1300, 1400, 1500, 1300, 1400, 1500],
            "Inhabitants": [5, 7, 9, 1, 2, 40],
        }
    )


@pytest.fixture
def yautepec_path() -> Path:
    path = DATA / "yautepec.csv"
    if not path.exists():
        pytest.skip("yautepec.csv is not available")
    return path


@pytest.fixture
def viabundus_path() -> Path:
    path = DATA / "viabundus.csv"
    if not path.exists():
        pytest.skip("viabundus.csv is not available")
    return path
