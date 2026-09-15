from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from modelcity import (
    as_settlements,
    kde,
    rank_size,
    read_yautepec,
    settlement_churn,
    settlement_duration,
    size_summary,
    size_totals,
    urban_share,
)

# Produced by R:
#   x <- c(120, 250, 260, 400, 900, 1500, 1520, 3000)
#   density(x, bw = 150, n = 9, from = 0, to = 3200)$y
R_DENSITY = np.array(
    [
        4.079674610e-04, 8.086220227e-04, 2.765801570e-04, 1.242610673e-04,
        5.543644664e-04, 3.281538809e-06, 1.124335153e-07, 1.367111848e-04,
        1.367275439e-04,
    ]
)


def test_kde_matches_r_density_with_the_same_bandwidth():
    """scipy scales the bandwidth by the sample spread; R does not.

    Remaining differences come from R binning onto a grid and using an FFT,
    where this evaluates the kernel directly.
    """
    result = kde([120, 250, 260, 400, 900, 1500, 1520, 3000], bw=150, n=9, lower=0, upper=3200)
    assert np.allclose(result["density"].to_numpy(), R_DENSITY, atol=1e-6)


def test_kde_without_an_explicit_range_extends_by_three_bandwidths():
    result = kde([100, 200], bw=10)
    assert result["year"].min() == pytest.approx(100 - 3 * 10)
    assert result["year"].max() == pytest.approx(200 + 3 * 10)


def test_kde_handles_a_single_observation():
    result = kde([500], bw=50)
    assert len(result) == 512
    assert np.isfinite(result["density"]).all()


def test_duration_spans_first_start_to_last_end(interval_frame):
    settlements = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    spans = settlement_duration(settlements).set_index("id")
    # Site A runs 3449-3050 and 3049-2450 BP, so 3449 - 2450 years.
    assert spans.loc["A", "duration"] == pytest.approx(999)
    assert spans.loc["B", "duration"] == pytest.approx(3049 - 2450)


def test_yautepec_median_duration_is_unchanged(yautepec_path):
    """Golden value: the R analysis reports a median site duration of 299 years."""
    settlements = read_yautepec(yautepec_path)
    spans = settlement_duration(settlements)
    assert spans["duration"].median() == pytest.approx(299.0)
    assert len(spans) == 782


def test_churn_is_the_sum_and_net_is_the_difference(yautepec_path):
    churn = settlement_churn(read_yautepec(yautepec_path), bw=150, n=64)
    assert np.allclose(churn["churn"], churn["foundations"] + churn["abandonments"])
    assert np.allclose(churn["net"], churn["foundations"] - churn["abandonments"])


def test_size_totals_counts_only_occupied_settlements():
    frame = pd.DataFrame(
        {
            "id": ["a", "b", "c"],
            "size": [10.0, 0.0, 5.0],
            "start": [100, 100, 100],
            "end": [200, 200, 200],
        }
    )
    settlements = as_settlements(frame, size_type="area", time_scale="CE")
    totals = size_totals(settlements)
    assert totals.loc[0, "total_size"] == pytest.approx(15.0)
    assert totals.loc[0, "n_settlements"] == 2


def test_size_summary_reports_the_quantile_band(interval_frame):
    settlements = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    summary = size_summary(settlements)
    assert (summary["q10"] <= summary["median"]).all()
    assert (summary["median"] <= summary["q90"]).all()
    assert (summary["max"] >= summary["q90"]).all()


def test_rank_size_ideal_is_the_zipf_line(interval_frame):
    settlements = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    ranked = rank_size(settlements)
    for _, period in ranked.groupby(["t_start", "t_end"]):
        largest = period["size"].max()
        assert np.allclose(period["ideal"], largest / period["rank"])
        assert period["rank"].tolist() == sorted(period["rank"].tolist())


def test_urban_share_defaults_follow_size_type(interval_frame, snapshot_frame):
    areas = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    assert "urban_20" in urban_share(areas).columns

    people = as_settlements(
        snapshot_frame, id="Nodes_ID", size="Inhabitants", start="Year",
        size_type="population", size_scale=1000,
    )
    assert "urban_1000" in urban_share(people).columns


def test_urban_share_percentages_are_consistent(interval_frame):
    settlements = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    shares = urban_share(settlements, thresholds=[20])
    expected = shares["urban_20"] / shares["total_size"] * 100
    assert np.allclose(shares["percent_20"], expected)


def test_custom_thresholds_are_honoured(interval_frame):
    settlements = as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )
    shares = urban_share(settlements, thresholds=[1, 2.5])
    assert {"urban_1", "urban_2.5"}.issubset(shares.columns)
