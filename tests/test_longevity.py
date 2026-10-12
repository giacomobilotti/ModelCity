from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.figure import Figure

from modelcity import (
    as_settlements,
    compare_survival,
    fit_persistence,
    fit_survival,
    impute_lifespan,
    mark_censoring,
    plot_longevity,
    sample_trajectories,
    summarize_persistence,
)
from modelcity.longevity import _prentice_gengamma
from modelcity.plots import ggplot_histogram_breaks, r_histogram_breaks
from modelcity.sampling import persistence_samples

# flexsurv dgengamma / pgengamma / qgengamma at mu=5.7, sigma=1.3, Q=0.9.
_MU, _SIGMA, _Q = 5.7, 1.3, 0.9
_GRID = np.array([1, 10, 50, 100, 400, 1000, 2500], dtype=float)
_PDF = np.array(
    [
        0.00738270842944,
        0.00481117787798,
        0.00299381746658,
        0.00217130986148,
        0.000698970231064,
        0.000160513423162,
        1.12667454036e-05,
    ]
)
_CDF = np.array(
    [
        0.0087306955358,
        0.0593612285093,
        0.206588328744,
        0.333384206787,
        0.70223193406,
        0.913098412736,
        0.992132313006,
    ]
)
_PROBABILITY = np.array([0.01, 0.1, 0.5, 0.9, 0.99])
_QUANTILE = np.array([1.17424315694, 19.2064875458, 195.977443125, 924.780140457, 2334.13954657])

# The same functions at Q=-0.4.
_Q_NEGATIVE = -0.4
_PDF_NEGATIVE = np.array(
    [
        1.8876158531e-09,
        0.000205239814757,
        0.00192534148983,
        0.00203381807929,
        0.000738796390786,
        0.000206522569792,
        4.08914342598e-05,
    ]
)
_CDF_NEGATIVE = np.array(
    [
        1.97809131682e-10,
        0.000518000733617,
        0.0503879237113,
        0.15335951214,
        0.534320333118,
        0.76954099735,
        0.908850042063,
    ]
)
_QUANTILE_NEGATIVE = np.array([24.5955436134, 74.3070277878, 356.483812564, 2300.14968495, 13888.5180848])

# A fixed censored sample and the flexsurvreg estimates fitted to it.
_TIMES = np.array([12, 40, 55, 80, 120, 200, 350, 500, 900, 1400, 30, 70, 150, 260, 600], dtype=float)
_EVENT = np.array([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1], dtype=float)
_FLEXSURV = {
    "exponential": {"rate": 0.0023075309},
    "weibull": {"shape": 0.7231755563, "scale": 410.1101465952},
    "lognormal": {"meanlog": 5.3268543321, "sdlog": 1.6156842728},
    "gompertz": {"shape": -0.0017685274, "rate": 0.0039003585},
    "gamma": {"shape": 0.6625906133, "rate": 0.0013829873},
    "llogis": {"shape": 1.0206327915, "scale": 200.6260404365},
    "gengamma": {"mu": 4.7414138448, "sigma": 1.5913044135, "Q": -0.8054162253},
}
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


_FLEXSURV_LOGLIK = {
    "exponential": -77.78734901,
    "weibull": -76.71458732,
    "lognormal": -75.86004676,
    "gompertz": -76.46455794,
    "gamma": -76.98202762,
    "llogis": -76.19374851,
    "gengamma": -75.64381162,
}


def test_censoring_follows_the_terminus_rule(phased):
    trajectories = sample_trajectories(phased, 30, seed=10)
    marked = mark_censoring(persistence_samples(trajectories, terminus=1751))
    reached = marked[marked["id"].isin(["A", "C"])]
    ended = marked[marked["id"] == "B"]
    assert (reached["event"] == 0).all()
    assert (ended["event"] == 1).all()
    assert reached["norm_pers"].eq(1).all()


def test_prentice_gengamma_matches_flexsurv():
    positive = _prentice_gengamma(_MU, _SIGMA, _Q)
    negative = _prentice_gengamma(_MU, _SIGMA, _Q_NEGATIVE)
    assert np.allclose(positive.pdf(_GRID), _PDF, rtol=1e-8, atol=1e-12)
    assert np.allclose(positive.cdf(_GRID), _CDF, rtol=1e-8, atol=1e-12)
    assert np.allclose(positive.ppf(_PROBABILITY), _QUANTILE, rtol=1e-8, atol=1e-8)
    assert np.allclose(negative.pdf(_GRID), _PDF_NEGATIVE, rtol=1e-8, atol=1e-14)
    assert np.allclose(negative.cdf(_GRID), _CDF_NEGATIVE, rtol=1e-8, atol=1e-12)
    assert np.allclose(negative.ppf(_PROBABILITY), _QUANTILE_NEGATIVE, rtol=1e-8, atol=1e-6)


def test_censored_fits_agree_with_flexsurv():
    fitted = fit_survival(_TIMES, _EVENT)
    assert set(fitted) == set(_FLEXSURV)
    for name, expected in _FLEXSURV.items():
        fit = fitted[name]
        assert list(fit.parameters) == list(expected)
        for parameter, value in expected.items():
            assert fit.parameters[parameter] == pytest.approx(value, rel=1e-3, abs=1e-3)
        assert fit.log_likelihood == pytest.approx(_FLEXSURV_LOGLIK[name], abs=1e-3)
        assert fit.n_params == len(expected)
        assert fit.aic == pytest.approx(-2 * _FLEXSURV_LOGLIK[name] + 2 * len(expected), abs=1e-2)


def test_comparison_ranks_by_aic_and_weights_sum_to_one():
    comparison = compare_survival(_TIMES, _EVENT)
    assert comparison["Rank"].tolist() == list(range(1, len(comparison) + 1))
    assert comparison["AIC"].is_monotonic_increasing
    assert comparison["Delta_AIC"].iloc[0] == pytest.approx(0.0)
    assert comparison["Weight"].sum() == pytest.approx(1.0)
    assert (comparison["Weight"] > 0).all()
    loglogistic = compare_survival(_TIMES, _EVENT, distributions=["loglogistic", "gengamma"])
    assert loglogistic["Model"].tolist()
    assert set(loglogistic["Model"]) == {"llogis", "gengamma"}


def test_imputation_is_a_reproducible_whole_year_at_least_as_long():
    observed = np.array([10.0, 20.0, 30.0])
    first = impute_lifespan(observed, 5.0, 1.2, 0.8, seed=1234)
    second = impute_lifespan(observed, 5.0, 1.2, 0.8, seed=1234)
    other_seed = impute_lifespan(observed, 5.0, 1.2, 0.8, seed=7)
    other_coefficients = impute_lifespan(observed, 6.5, 0.4, -0.3, seed=1234)

    assert np.array_equal(first, second)
    assert not np.array_equal(first, other_seed)
    assert not np.array_equal(first, other_coefficients)
    assert np.all(first == np.floor(first))
    assert np.all(first >= observed)


def test_non_finite_conditional_quantile_is_missing():
    imputed = impute_lifespan([1e12], _MU, _SIGMA, _Q, seed=1)
    assert np.isnan(imputed[0])


def test_non_positive_lifespans_are_rejected():
    with pytest.raises(ValueError, match="strictly positive"):
        fit_survival([10, 0, 5], [1, 1, 1], distributions="exponential")
    with pytest.raises(ValueError, match="strictly positive"):
        impute_lifespan([10, -1], 5.0, 1.0, 0.5)
    with pytest.raises(ValueError, match="unknown distribution"):
        fit_survival(_TIMES, _EVENT, distributions="pareto")


def test_persistence_fits_one_row_per_converged_iteration(phased):
    trajectories = sample_trajectories(phased, 4, seed=3, terminus=1751)
    samples = persistence_samples(trajectories, terminus=1751)
    fitted = fit_persistence(samples, distribution="gengamma")
    assert fitted["iteration"].tolist() == [1, 2, 3, 4]
    assert {"mu", "sigma", "Q", "aic"}.issubset(fitted.columns)
    assert np.isfinite(fitted[["mu", "sigma", "Q"]]).all().all()
    assert (fitted["sigma"] > 0).all()


def test_summary_keeps_the_share_from_before_imputation():
    samples = pd.DataFrame(
        {
            "id": ["A", "A", "B"],
            "iteration": [1, 2, 1],
            "available_pers": [100.0, 100.0, 80.0],
            "persistence": [250.0, 50.0, 80.0],
            "norm_pers": [1.0, 0.5, 1.0],
            "start": [300.0, 280.0, 200.0],
            "end": [100.0, 150.0, 120.0],
        }
    )
    summary = summarize_persistence(samples).set_index("id")
    assert summary.loc["A", "mean"] == pytest.approx(150.0)
    assert summary.loc["A", "mean_n"] == pytest.approx(0.75)
    assert summary.loc["B", "mean_n"] == pytest.approx(1.0)
    assert list(summary.columns) == [
        "mean",
        "median",
        "q025",
        "q975",
        "q25",
        "q75",
        "mean_n",
        "median_n",
        "q025_n",
        "q975_n",
        "q25_n",
        "q75_n",
        "start",
        "end",
    ]


def test_summary_skips_a_missing_share():
    samples = pd.DataFrame(
        {
            "id": ["A", "A"],
            "persistence": [10.0, 40.0],
            "norm_pers": [np.nan, 0.4],
            "start": [1.0, 1.0],
            "end": [1.0, 1.0],
        }
    )
    summary = summarize_persistence(samples).iloc[0]
    assert summary["mean"] == pytest.approx(25.0)
    assert summary["mean_n"] == pytest.approx(0.4)


def test_histogram_breaks_match_ggplot_and_base_r():
    durations = np.array([1, 2, 5, 10, 20, 50, 100, 400, 1200], dtype=float)
    breaks = ggplot_histogram_breaks(durations, bins=40)
    assert breaks[0] == pytest.approx(-14.371794871794872)
    assert breaks[-1] == pytest.approx(1215.3717948717949)
    assert len(breaks) == 41
    short = np.array([12, 40, 55, 80, 120], dtype=float)
    assert ggplot_histogram_breaks(short, bins=5) == pytest.approx([-1.5, 25.5, 52.5, 79.5, 106.5, 133.5])
    integers = np.arange(40, dtype=float)
    assert ggplot_histogram_breaks(integers, bins=40)[0] == pytest.approx(-0.5)
    pretty = r_histogram_breaks(np.array([3.2, 2480.7]), bins=40)
    assert pretty[0] == pytest.approx(0)
    assert pretty[-1] == pytest.approx(2500)
    assert np.diff(pretty) == pytest.approx(50)
    assert r_histogram_breaks(np.array([0.0, 0.0]), bins=40) == pytest.approx([-1.0, 0.0])


def test_longevity_plot_draws_the_fitted_curve():
    samples = pd.DataFrame({"persistence": [10.0, 20.0, 40.0, 80.0, 160.0]})
    figure = plot_longevity(samples, coefficients={"mu": 4.0, "sigma": 1.0, "Q": 0.8})
    assert isinstance(figure, Figure)
    assert len(figure.axes[0].lines) == 1
    plt.close(figure)
    with pytest.raises(ValueError, match="mu"):
        plot_longevity(samples, coefficients={"mu": 4.0})
