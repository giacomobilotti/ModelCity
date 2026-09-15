from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pytest
from matplotlib.figure import Figure

from modelcity import (
    as_settlements,
    plot_churn,
    plot_duration,
    plot_events,
    plot_rank_size,
    plot_size_summary,
    plot_size_totals,
    plot_urban_share,
    read_viabundus,
    read_yautepec,
    year_display,
)

PLOTS = [
    plot_duration,
    plot_events,
    plot_churn,
    plot_size_totals,
    plot_size_summary,
    plot_rank_size,
    plot_urban_share,
]


@pytest.fixture
def areas(interval_frame):
    return as_settlements(
        interval_frame, id="Sitio", size="Area", size_type="area", time_scale="BP"
    )


@pytest.fixture
def people(snapshot_frame):
    return as_settlements(
        snapshot_frame,
        id="Nodes_ID",
        size="Inhabitants",
        start="Year",
        name="Name",
        size_type="population",
        size_scale=1000,
    )


@pytest.mark.parametrize("plot", PLOTS, ids=lambda fn: fn.__name__)
def test_every_plot_works_for_area_data(plot, areas):
    figure = plot(areas)
    assert isinstance(figure, Figure)
    plt.close(figure)


@pytest.mark.parametrize("plot", PLOTS, ids=lambda fn: fn.__name__)
def test_every_plot_works_for_population_data(plot, people):
    """The same calls must work whichever size measure the dataset uses."""
    figure = plot(people)
    assert isinstance(figure, Figure)
    plt.close(figure)


def test_axis_labels_follow_the_size_type(areas, people):
    area_figure = plot_size_totals(areas)
    assert "area" in area_figure.axes[0].get_ylabel().lower()
    assert "[ha]" in area_figure.axes[0].get_ylabel()
    plt.close(area_figure)

    people_figure = plot_size_totals(people)
    assert "population" in people_figure.axes[0].get_ylabel().lower()
    assert "[inhabitants]" in people_figure.axes[0].get_ylabel()
    plt.close(people_figure)


def test_size_totals_has_a_second_axis_for_counts(areas):
    figure = plot_size_totals(areas)
    assert len(figure.axes) == 2
    assert "Number of" in figure.axes[1].get_ylabel()
    plt.close(figure)


def test_duration_title_reports_the_median(areas):
    figure = plot_duration(areas)
    assert "median" in figure.axes[0].get_title()
    plt.close(figure)


def test_year_axis_is_labelled_in_bp_by_default(areas):
    figure = plot_events(areas)
    ax = figure.axes[0]
    assert ax.get_xlabel() == "Year BP"

    # Ticks land on round BP values rather than on round astronomical ones,
    # which would label the axis 3450, 2950, 2450.
    figure.canvas.draw()
    labels = [float(text.get_text()) for text in ax.get_xticklabels() if text.get_text()]
    assert len(labels) > 2
    step = labels[0] - labels[1]
    assert step % 50 == 0
    assert all(value % step == 0 for value in labels)
    # Time runs left to right, so BP counts down.
    assert labels == sorted(labels, reverse=True)
    plt.close(figure)


def test_year_axis_can_be_switched_back_to_bc_ad(areas):
    with year_display("BCAD"):
        figure = plot_events(areas)
    assert figure.axes[0].get_xlabel() == "Year (BC/AD)"
    plt.close(figure)


def test_duration_axis_stays_in_elapsed_years(areas):
    """Durations are spans, not dates, so no display scale applies to them."""
    figure = plot_duration(areas)
    assert figure.axes[0].get_xlabel() == "Duration (years)"
    plt.close(figure)


def test_display_scale_cannot_change_the_numbers(areas):
    """Switching the labels must move no data, only how it is written."""
    baseline = plot_size_totals(areas).axes[0].lines
    baseline_points = [line.get_xydata().tolist() for line in baseline]
    with year_display("BCAD"):
        switched = plot_size_totals(areas).axes[0].lines
        switched_points = [line.get_xydata().tolist() for line in switched]
    assert baseline_points == switched_points
    plt.close("all")


def test_plots_run_on_the_real_datasets(yautepec_path, viabundus_path):
    for settlements in (read_yautepec(yautepec_path), read_viabundus(viabundus_path)):
        for plot in PLOTS:
            figure = plot(settlements)
            assert isinstance(figure, Figure)
            plt.close(figure)
