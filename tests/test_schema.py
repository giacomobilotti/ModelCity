from __future__ import annotations

import pandas as pd
import pytest

from modelcity import SchemaError, as_settlements
from modelcity.schema import resolve_roles, suggest_column


def test_canonical_names_resolve_without_being_named(interval_frame):
    """Columns already called `start` and `end` need no mapping."""
    roles = resolve_roles(
        interval_frame.columns, id="Sitio", size="Area", size_type="area"
    )
    assert roles.start == "start"
    assert roles.end == "end"


def test_case_only_difference_is_accepted():
    frame = pd.DataFrame({"ID": [1], "Size": [2.0], "Start": [100], "End": [50]})
    roles = resolve_roles(frame.columns, size_type="area")
    assert (roles.id, roles.size, roles.start, roles.end) == ("ID", "Size", "Start", "End")


def test_missing_size_role_names_the_role_and_suggests_a_column(interval_frame):
    with pytest.raises(SchemaError) as excinfo:
        as_settlements(interval_frame, id="Sitio", size_type="area", time_scale="BP")

    message = str(excinfo.value)
    assert "role 'size'" in message
    assert "did you mean 'Area'?" in message
    assert "available:" in message


def test_error_lists_every_unmapped_role_at_once(interval_frame):
    """One pass should surface every problem, not just the first."""
    with pytest.raises(SchemaError) as excinfo:
        as_settlements(interval_frame, size_type="area", time_scale="BP")

    message = str(excinfo.value)
    assert "role 'id'" in message
    assert "role 'size'" in message


def test_fix_line_is_copy_pasteable(interval_frame):
    with pytest.raises(SchemaError) as excinfo:
        as_settlements(interval_frame, size_type="area", time_scale="BP")

    fix = [line for line in str(excinfo.value).splitlines() if line.strip().startswith("fix:")]
    assert len(fix) == 1
    # Roles that resolved keep their real column; the failures use suggestions.
    assert 'size="Area"' in fix[0]
    assert 'start="start"' in fix[0]
    assert 'size_type="area"' in fix[0]


def test_supplied_column_that_does_not_exist_is_reported(interval_frame):
    with pytest.raises(SchemaError) as excinfo:
        as_settlements(
            interval_frame, id="Sitio", size="Hectares", size_type="area", time_scale="BP"
        )

    message = str(excinfo.value)
    assert "'Hectares' was supplied for role 'size' but is not in the data" in message
    assert "did you mean 'Area'?" in message


def test_size_type_is_required_and_checked(interval_frame):
    with pytest.raises(SchemaError, match="size_type is required"):
        as_settlements(interval_frame, id="Sitio", size="Area", time_scale="BP")

    with pytest.raises(SchemaError, match="must be one of"):
        as_settlements(
            interval_frame, id="Sitio", size="Area", size_type="hectares", time_scale="BP"
        )


def test_time_scale_is_checked(interval_frame):
    with pytest.raises(SchemaError, match="time_scale must be"):
        as_settlements(
            interval_frame, id="Sitio", size="Area", size_type="area", time_scale="bp"
        )


@pytest.mark.parametrize(
    "role, columns, expected",
    [
        ("size", ["Sitio", "Area", "start"], "Area"),
        ("size", ["Nodes_ID", "Inhabitants"], "Inhabitants"),
        ("id", ["Sitio", "Area"], "Sitio"),
        ("id", ["Nodes_ID", "Year"], "Nodes_ID"),
        ("group", ["region", "Year"], "region"),
        ("size", ["colour", "shape"], None),
    ],
)
def test_suggestions_use_domain_aliases(role, columns, expected):
    """`size` and `Area` share no characters, so plain fuzzy matching is not enough."""
    assert suggest_column(role, columns) == expected


def test_non_dataframe_input_is_rejected():
    with pytest.raises(SchemaError, match="must be a pandas DataFrame"):
        as_settlements([1, 2, 3], id="a", size="b", start="c", size_type="area")
