from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from modelcity import (
    SchemaError,
    as_settlements,
    period_summary,
    persistence_samples,
    persistence_summary,
    plot_period_growth,
    plot_persistence,
    plot_site_growth,
    plot_trajectories,
    sample_trajectories,
    site_summary,
    to_bp,
)
from run_unified import read_yautepec

METHODS = [
    dict(start_method="unif", end_method="unif", peak_method="unif"),
    dict(start_method="exp", end_method="exp", peak_method="norm"),
    dict(start_method="exp", end_method="unif", peak_method="norm"),
]


@pytest.fixture
def phased():
    """A spans three contiguous phases, B one early phase, C the last phase."""
    frame = pd.DataFrame(
        {
            "Sitio": ["A", "A", "A", "B", "C"],
            "Area": [1.0, 4.0, 10.0, 2.0, 3.0],
            "start": [3049, 2449, 2049, 3049, 2049],
            "end": [2450, 2050, 1751, 2450, 1751],
        }
    )
    return as_settlements(frame, id="Sitio", size="Area", size_type="area", time_scale="BP")


def _bp(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.assign(
        year_bp=to_bp(frame["year"]), start_bp=to_bp(frame["t_start"]), end_bp=to_bp(frame["t_end"])
    )


@pytest.mark.parametrize("methods", METHODS)
def test_dates_are_whole_years_inside_their_phase(phased, methods):
    trajectories = _bp(sample_trajectories(phased, 400, seed=1, **methods))
    assert (trajectories["year_bp"] % 1 == 0).all()
    assert (trajectories["year_bp"] <= trajectories["start_bp"]).all()
    assert (trajectories["year_bp"] >= trajectories["end_bp"]).all()


@pytest.mark.parametrize("methods", METHODS[:2])
def test_multi_phase_draws_run_forward_in_time(phased, methods):
    trajectories = sample_trajectories(phased, 400, seed=1, **methods)
    a = trajectories[trajectories["id"] == "A"]
    for _, draw in a.groupby("iteration"):
        assert np.all(np.diff(draw["year"].to_numpy()) >= 0)


def test_steps_pair_each_date_with_the_recorded_size(phased):
    # Every end in A's last phase (2049-1751 BP) falls at or after 2049 BP.
    trajectories = sample_trajectories(phased, 3, seed=1, terminus=2049)
    site = trajectories[(trajectories["id"] == "A") & (trajectories["iteration"] == 1)]
    assert site["event"].tolist() == ["foundation", "peak", "peak", "peak", "abandonment"]
    assert site["size"].tolist() == [0.1, 1.0, 4.0, 10.0, 10.0]
    assert site["phase"].tolist() == [0, 0, 1, 2, 2]


def test_terminus_sets_abandonment_to_zero_per_draw(phased):
    # R's default terminus of 1 BP zeroes every ending.
    default = sample_trajectories(phased, 50, seed=2)
    assert (default.loc[default["event"] == "abandonment", "size"] == 0).all()

    # With the terminus inside C's phase, only draws ending before it are zeroed.
    trajectories = _bp(sample_trajectories(phased, 2000, seed=2, terminus=1900))
    ends = trajectories[(trajectories["event"] == "abandonment") & (trajectories["id"] == "C")]
    assert ((ends["size"] == 0) == (ends["year_bp"] > 1900)).all()
    assert 0 < (ends["size"] == 0).mean() < 1


def test_seed_makes_draws_reproducible(phased):
    first = sample_trajectories(phased, 20, seed=7)
    second = sample_trajectories(phased, 20, seed=7)
    pd.testing.assert_frame_equal(first, second)


def test_exponential_start_matches_rs_floored_truncated_exponential():
    frame = pd.DataFrame({"id": ["a"], "size": [1.0], "start": [1000], "end": [0]})
    settlements = as_settlements(frame, size_type="area", time_scale="BP")
    trajectories = _bp(sample_trajectories(settlements, 40000, start_method="exp", perc=5, seed=3))
    delay = 1000 - trajectories.loc[trajectories["event"] == "foundation", "year_bp"]
    assert (delay % 1 == 0).all() and delay.min() >= 0 and delay.max() < 1000
    # floor() of an exponential with mean 200, truncated at 1000.
    k = np.arange(1000)
    weights = np.exp(-k / 200.0) * (1 - np.exp(-1 / 200.0))
    assert delay.mean() == pytest.approx((k * weights).sum() / weights.sum(), rel=0.02)


def test_single_phase_exponential_end_follows_the_foundation():
    """R draws the single-phase end as a short delay after the foundation."""
    frame = pd.DataFrame({"id": ["a"], "size": [1.0], "start": [1000], "end": [0]})
    settlements = as_settlements(frame, size_type="area", time_scale="BP")
    trajectories = _bp(sample_trajectories(settlements, 20000, start_method="exp", end_method="exp", seed=4))
    first = trajectories.loc[trajectories["event"] == "foundation", "year_bp"].to_numpy()
    last = trajectories.loc[trajectories["event"] == "abandonment", "year_bp"].to_numpy()
    assert np.median(first - last) < 0.25 * np.median(last)


def test_single_phase_normal_peak_can_follow_the_abandonment_and_swaps_growth():
    frame = pd.DataFrame({"id": ["a"], "size": [5.0], "start": [1000], "end": [0]})
    settlements = as_settlements(frame, size_type="area", time_scale="BP")
    trajectories = sample_trajectories(settlements, 4000, peak_method="norm", terminus=1, seed=5)
    peak = trajectories[trajectories["event"] == "peak"].reset_index(drop=True)
    end = trajectories[trajectories["event"] == "abandonment"].reset_index(drop=True)
    late = peak["year"] > end["year"]
    assert late.any()
    # Sorted by date the abandonment comes second, so the peak's row carries the
    # rate into the abandonment (from 0.1 to zero, -1) as R writes it back.
    assert np.allclose(peak.loc[late, "growth_cagr"], -1.0)


def test_growth_matches_r_formula():
    frame = pd.DataFrame({"id": ["a", "a"], "size": [1.0, 4.0], "start": [3000, 2000], "end": [2001, 1000]})
    settlements = as_settlements(frame, size_type="area", time_scale="BP")
    trajectories = sample_trajectories(settlements, 50, terminus=1, seed=6)
    for _, draw in trajectories.groupby("iteration"):
        years, sizes = draw["year"].to_numpy(), draw["size"].to_numpy()
        interval = np.maximum(np.diff(years), 1.0)
        expected = (sizes[1:] / sizes[:-1]) ** (1 / interval) - 1
        assert np.isnan(draw["growth_cagr"].iloc[0])
        assert np.allclose(draw["growth_cagr"].iloc[1:], expected)
        assert np.allclose(draw["growth_linear"].iloc[1:], np.diff(sizes) / interval)


def test_site_summary_groups_by_phase_and_keeps_rs_missing_medians(phased):
    trajectories = sample_trajectories(phased, 30, seed=8, terminus=1751)
    summary = site_summary(trajectories)
    assert len(summary) == 3 + 1 + 1
    a = summary[summary["id"] == "A"]
    assert a["t_start"].is_monotonic_decreasing
    # The foundation's missing rate blanks the earliest phase, as median() does in R.
    assert np.isnan(a["median"].iloc[-1]) and np.isnan(a["mean"].iloc[-1])
    assert a["median"].iloc[:-1].notna().all()
    assert np.isfinite(a["q975"].iloc[-1])
    assert (summary["q025_year"] <= summary["q975_year"]).all()


def test_period_summary_skips_missing_values(phased):
    trajectories = sample_trajectories(phased, 30, seed=9)
    periods = period_summary(trajectories)
    assert len(periods) == phased.n_periods
    assert periods["median"].notna().all()
    assert periods["t_start"].is_monotonic_decreasing


def test_persistence_follows_rs_terminus_rule(phased):
    trajectories = sample_trajectories(phased, 200, seed=10)
    samples = persistence_samples(trajectories, terminus=1751)
    reached = samples[samples["id"].isin(["A", "C"])]
    assert np.allclose(reached["norm_pers"], 1.0)
    assert (reached["end"] == 1751).all()
    assert np.allclose(reached["persistence"], reached["available_pers"])

    lost = samples[samples["id"] == "B"]
    assert (lost["norm_pers"] < 1).all()
    assert (lost["persistence"] >= 1).all()
    assert np.allclose(lost["available_pers"], lost["start"] - 1751)

    summary = persistence_summary(trajectories, terminus=1751).set_index("id")
    assert summary.loc["B", "q025"] <= summary.loc["B", "median"] <= summary.loc["B", "q975"]


def test_persistence_needs_a_time_scale_once_metadata_is_lost(phased):
    trajectories = sample_trajectories(phased, 5, seed=10)
    stripped = trajectories.copy()
    stripped.attrs = {}
    with pytest.raises(ValueError, match="time_scale"):
        persistence_samples(stripped, terminus=1751)
    pd.testing.assert_frame_equal(
        persistence_samples(stripped, terminus=1751, time_scale="BP"),
        persistence_samples(trajectories, terminus=1751),
    )


def test_invalid_inputs_are_rejected(phased, snapshot_frame):
    with pytest.raises(ValueError, match="start_method"):
        sample_trajectories(phased, 5, start_method="gamma")
    with pytest.raises(ValueError, match="perc"):
        sample_trajectories(phased, 5, perc=0)

    duplicated = pd.DataFrame({"id": ["a", "a"], "size": [1.0, 2.0], "start": [0, 0], "end": [100, 50]})
    with pytest.raises(SchemaError, match="uplicates"):
        sample_trajectories(as_settlements(duplicated, size_type="area"), 5)

    instant = pd.DataFrame({"id": ["a"], "size": [1.0], "start": [100], "end": [100]})
    with pytest.raises(SchemaError, match="older than end"):
        sample_trajectories(as_settlements(instant, size_type="area"), 5)

    people = as_settlements(
        snapshot_frame, id="Nodes_ID", size="Inhabitants", start="Year", size_type="population"
    )
    with pytest.raises(SchemaError, match="snapshot"):
        sample_trajectories(people, 5)


def test_yautepec_sampling_covers_every_site(yautepec_path):
    settlements = read_yautepec(yautepec_path)
    trajectories = sample_trajectories(settlements, 5, terminus=510, seed=11)
    assert trajectories["id"].nunique() == settlements.n_settlements
    assert len(trajectories) == 5 * (len(settlements.data) + 2 * settlements.n_settlements)


@pytest.mark.parametrize(
    "draw",
    [
        lambda t: plot_site_growth(t, "A"),
        lambda t: plot_trajectories(t, "A"),
        lambda t: plot_trajectories(t, "A", value="area", seed=1),
        lambda t: plot_trajectories(t, "A", value="growth_linear"),
        lambda t: plot_period_growth(t),
        lambda t: plot_persistence(t, terminus=1751, normalised=True),
    ],
)
def test_sampling_plots_return_figures(phased, draw):
    trajectories = sample_trajectories(phased, 40, seed=12, terminus=1751)
    figure = draw(trajectories)
    assert isinstance(figure, Figure)
    plt.close(figure)


def test_site_growth_plot_needs_more_than_one_phase(phased):
    trajectories = sample_trajectories(phased, 5, seed=13)
    with pytest.raises(ValueError, match="single phase"):
        plot_site_growth(trajectories, "B")
    with pytest.raises(ValueError, match="no trajectories"):
        plot_site_growth(trajectories, "missing")
