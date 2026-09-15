from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from modelcity import SchemaError, as_settlements, year_display


def build_interval(frame, **kwargs):
    options = dict(id="Sitio", size="Area", size_type="area", time_scale="BP")
    options.update(kwargs)
    return as_settlements(frame, **options)


def test_bp_input_becomes_increasing_calendar_years(interval_frame):
    settlements = build_interval(interval_frame)
    assert (settlements.data["t_start"] <= settlements.data["t_end"]).all()
    # 3449 BP is the earliest start in the fixture.
    assert settlements.data["t_start"].min() == pytest.approx(1950 - 3449)


def test_interval_data_is_detected(interval_frame):
    settlements = build_interval(interval_frame)
    assert settlements.time.record_type == "interval"
    assert not settlements.is_snapshot


def test_snapshot_data_is_detected_and_binned(snapshot_frame):
    settlements = as_settlements(
        snapshot_frame,
        id="Nodes_ID",
        size="Inhabitants",
        start="Year",
        name="Name",
        size_type="population",
        size_scale=1000,
    )
    assert settlements.is_snapshot

    first = settlements.data[settlements.data["t_observed"] == 1400].iloc[0]
    assert (first["t_start"], first["t_end"]) == (1350.0, 1450.0)


def test_size_scale_brings_values_into_real_units(snapshot_frame):
    settlements = as_settlements(
        snapshot_frame,
        id="Nodes_ID",
        size="Inhabitants",
        start="Year",
        size_type="population",
        size_scale=1000,
    )
    assert settlements.data["size"].max() == 40_000


def test_size_type_drives_units_and_thresholds(interval_frame, snapshot_frame):
    areas = build_interval(interval_frame)
    assert areas.size.unit == "ha"
    assert areas.size.default_thresholds == (20.0, 40.0)
    assert "ha" in areas.size_label

    people = as_settlements(
        snapshot_frame, id="Nodes_ID", size="Inhabitants", start="Year",
        size_type="population",
    )
    assert people.size.unit == "inhabitants"
    assert people.size.default_thresholds == (1000.0, 2000.0)


def test_reversed_periods_are_reported_with_an_explanation(interval_frame):
    """Swapping start and end on BP data is a common mistake worth naming."""
    swapped = interval_frame.rename(columns={"start": "end", "end": "start"})
    with pytest.raises(SchemaError) as excinfo:
        build_interval(swapped)
    assert "end before they start" in str(excinfo.value)
    assert "BP counts backwards" in str(excinfo.value)


def test_negative_sizes_are_rejected(interval_frame):
    frame = interval_frame.copy()
    frame.loc[0, "Area"] = -5.0
    with pytest.raises(SchemaError, match="negative values"):
        build_interval(frame)


def test_unreadable_size_column_is_reported(interval_frame):
    frame = interval_frame.copy()
    frame["Area"] = ["small", "large", "huge", "tiny"]
    with pytest.raises(SchemaError, match="could not be read as numbers"):
        build_interval(frame)


def test_extra_columns_are_preserved(interval_frame):
    frame = interval_frame.copy()
    frame["excavated"] = [True, False, True, False]
    settlements = build_interval(frame)
    assert "excavated" in settlements.data.columns


def test_name_falls_back_to_the_identifier(interval_frame):
    settlements = build_interval(interval_frame)
    assert (settlements.data["name"] == settlements.data["id"]).all()


def test_periods_are_ordered_and_labelled(interval_frame):
    settlements = build_interval(interval_frame)
    periods = settlements.periods()
    assert periods["t_start"].is_monotonic_increasing
    # Labels default to BP, restating the years the source file was written in.
    assert periods["label"].iloc[0] == "3449-3050 BP"
    with year_display("BCAD"):
        assert settlements.periods()["label"].iloc[0] == "1500-1101 BC"


def test_years_bp_restates_the_source_years(interval_frame):
    settlements = build_interval(interval_frame)
    restated = settlements.years_bp()
    first = restated.iloc[0]
    assert (first["start_bp"], first["end_bp"]) == (3449, 3050)
    # BP counts backwards, so the start is the larger number of the pair.
    assert (restated["start_bp"] > restated["end_bp"]).all()


def test_repr_summarises_the_dataset(interval_frame):
    text = repr(build_interval(interval_frame))
    assert "3 settlements" in text
    assert "area in ha" in text
