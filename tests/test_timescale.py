from __future__ import annotations

import numpy as np
import pytest

from modelcity.timescale import (
    bcad_to_bp,
    bp_to_ce,
    ce_to_bp,
    format_period,
    format_year,
    get_year_display,
    set_year_display,
    snapshot_bounds,
    to_astronomical,
    to_bp,
    year_display,
)


# -- parity with the R implementation ---------------------------------------
#
# Expected values were produced by running the project's R `bcad_to_bp()` over
# these vectors, so a divergence here means the port has drifted.

R_BCAD_TO_BP = [
    (-5000, 6949), (-1499, 3448), (-100, 2049), (-2, 1951), (-1, 1950),
    (1, 1949), (2, 1948), (100, 1850), (1500, 450), (1649, 301),
    (1949, 1), (1950, 0), (2000, 50),
]

R_BP_TO_BCAD = [
    (0, 1950), (1, 1949), (300, 1650), (1949, 1), (1950, -1),
    (1951, -2), (2500, -551), (3449, -1500), (5000, -3051),
]


@pytest.mark.parametrize("bcad, expected_bp", R_BCAD_TO_BP)
def test_bcad_to_bp_matches_r(bcad, expected_bp):
    assert bcad_to_bp(bcad) == expected_bp


@pytest.mark.parametrize("bp, expected_bcad", R_BP_TO_BCAD)
def test_bp_to_bcad_matches_r(bp, expected_bcad):
    assert bcad_to_bp(bp, bc_to_bp=False) == expected_bcad


def test_conversion_runs_over_vectors_and_keeps_missing_values():
    """`NaN` passes through as `NA` does in R, rather than converting."""
    result = bcad_to_bp(np.array([-1499, np.nan, 1500]))
    assert result[0] == 3448
    assert np.isnan(result[1])
    assert result[2] == 450


def test_year_zero_is_rejected_like_in_r():
    with pytest.raises(ValueError, match="0 BC/AD is not a valid year"):
        bcad_to_bp(0)
    with pytest.raises(ValueError, match="0 BC/AD is not a valid year"):
        bcad_to_bp([-100, 0, 100])


def test_post_bomb_dates_are_rejected_like_in_r():
    with pytest.raises(ValueError, match="Post-bomb dates"):
        bcad_to_bp(-5, bc_to_bp=False)


def test_named_wrappers_delegate_to_the_r_port():
    assert ce_to_bp(-1499) == bcad_to_bp(-1499)
    assert bp_to_ce(3449) == bcad_to_bp(3449, bc_to_bp=False)


def test_bp_to_ce_skips_the_year_that_never_existed():
    assert bp_to_ce(0) == 1950
    assert bp_to_ce(1000) == 950
    assert bp_to_ce(1950) == -1
    assert bp_to_ce(1951) == -2


def test_bp_to_ce_round_trips():
    values = np.array([0, 500, 1949, 1950, 2500, 3449])
    assert np.allclose(ce_to_bp(bp_to_ce(values)), values)


# -- the internal astronomical scale ----------------------------------------

def test_astronomical_years_preserve_elapsed_time():
    """Differences must equal the difference of the BP values.

    An interval spanning the BC/AD boundary would come out a year too long if
    durations were computed on historical years, because of the missing zero.
    """
    start_bp, end_bp = 2049, 1751
    span = to_astronomical(end_bp) - to_astronomical(start_bp)
    assert span == start_bp - end_bp == 298

    # The historical numbering used for labels does differ, by exactly one year.
    label_span = bp_to_ce(end_bp) - bp_to_ce(start_bp)
    assert label_span == 299


def test_to_bp_inverts_to_astronomical_without_a_boundary_case():
    values = np.array([0, 431, 1950, 3449, 5000])
    assert np.allclose(to_bp(to_astronomical(values)), values)
    # Fractional midpoints survive, which the no-year-zero rule cannot promise.
    assert to_bp(to_astronomical(3449.5)) == 3449.5


@pytest.mark.parametrize("bp", [431, 1751, 1949, 1950, 2049, 3449])
def test_both_routes_to_bp_agree(bp):
    """Going through astronomical years must land where the R route lands.

    The internal scale reaches BP by plain subtraction while the R conversion
    steps over the missing year zero; they have to agree, or stored years and
    published labels would describe different dates.
    """
    assert to_bp(to_astronomical(bp)) == ce_to_bp(bp_to_ce(bp)) == bp


# -- display scale ----------------------------------------------------------

def test_bp_is_the_default_display():
    assert get_year_display() == "BP"
    assert format_year(-1499) == "3449 BP"
    assert format_year(1500) == "450 BP"


def test_display_can_be_switched_to_bc_ad():
    set_year_display("BCAD")
    assert format_year(1500) == "1500 AD"
    assert format_year(0) == "1 BC"
    assert format_year(-1499) == "1500 BC"
    assert format_year(-1499, suffix=False) == "-1500"


def test_display_can_be_overridden_per_call():
    assert format_year(-1499, display="BCAD") == "1500 BC"
    assert get_year_display() == "BP"


def test_year_display_context_manager_restores_the_previous_scale():
    with year_display("BCAD"):
        assert format_year(-1499) == "1500 BC"
    assert format_year(-1499) == "3449 BP"


def test_unknown_display_is_rejected():
    with pytest.raises(ValueError, match="year display must be"):
        set_year_display("BCE")


def test_period_labels_write_the_era_once():
    assert format_period(-1499, -1100) == "3449-3050 BP"
    with year_display("BCAD"):
        assert format_period(-1499, -1100) == "1500-1101 BC"
        # An era change has to be spelled out on both ends to stay readable.
        assert format_period(-49, 50) == "50 BC to 50 AD"


def test_format_year_ignores_missing_values():
    assert format_year(float("nan")) == ""


# -- snapshot binning -------------------------------------------------------

def test_snapshot_bounds_are_midpoints_between_observations():
    bounds = snapshot_bounds([1300, 1400, 1500, 1550])
    assert bounds[1400] == (1350.0, 1450.0)
    assert bounds[1500] == (1450.0, 1525.0)


def test_snapshot_bounds_extend_the_outer_observations():
    """The first and last bins mirror their single neighbouring gap."""
    bounds = snapshot_bounds([1300, 1400, 1500])
    assert bounds[1300] == (1250.0, 1350.0)
    assert bounds[1500] == (1450.0, 1550.0)


def test_snapshot_bounds_handle_degenerate_input():
    assert snapshot_bounds([]) == {}
    assert snapshot_bounds([1400]) == {1400.0: (1400.0, 1400.0)}
