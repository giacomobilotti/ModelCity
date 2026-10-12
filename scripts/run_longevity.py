"""Longevity of the example datasets, after chronological sampling.

Draws foundation and abandonment dates, fits a generalised gamma to each
iteration's censored lifespans, and imputes a total lifespan where the
settlement was still occupied at the persistence terminus.  Figures and tables
go to ``figures/<dataset>/longevity-python/``.

The area survey follows ``scripts/r/run-longevity-plots.R``: every site, 1000
iterations, sampling to 510 BP and persistence to 431 BP.  Draws come from
NumPy, so individual imputations differ from the R run.  The fitted
coefficients, the histogram bins, and the normalised share (taken before
imputation) follow the R script.  Large / Failed Large / Small labels are that
survey's rule (10 ha, and a phase starting at 510 BP or later) and are not
part of the library.

The population survey uses the same phase reading as ``scripts/run_sampling.py``:
each census is the period inferred around it, zero-population rows are dropped,
sampling stops at 1625 CE and persistence is measured to 1675 CE.  It has no
size-class figures.

    python scripts/run_longevity.py
    python scripts/run_longevity.py --dataset viabundus
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import modelcity as mc  # noqa: E402
from modelcity.longevity import _prentice_gengamma
from modelcity.plots import ggplot_histogram_breaks, r_histogram_breaks  # noqa: E402
from run_sampling import census_phases  # noqa: E402
from run_unified import DATA, FIGURES, ROOT, read_viabundus, read_yautepec  # noqa: E402

ITERATIONS = 1000
SEED = 42
SIZE_COLOURS = {"Large": "#CD0000", "Failed Large": "#6CA6CD", "Small": "#008B00"}
SIZE_ORDER = ("Large", "Failed Large", "Small")


def size_labels(path) -> pd.Series:
    """Size classes, keyed by site id.  Phase dates in the CSV are years BP."""
    sites = pd.read_csv(path)
    sites["Sitio"] = sites["Sitio"].astype(str)
    large = set(sites.loc[(sites["Area"] >= 10) & (sites["start"] <= 510), "Sitio"])
    failed = set(sites.loc[sites["Area"] >= 10, "Sitio"]) - large

    def label(identifier: str) -> str:
        if identifier in large:
            return "Large"
        if identifier in failed:
            return "Failed Large"
        return "Small"

    return label


def save(figure, destination, name: str) -> None:
    path = mc.save_figure(figure, destination / "{}.png".format(name))
    plt.close(figure)
    print("  wrote {}".format(Path(path).relative_to(ROOT)))


def plot_against_foundation(
    summary: pd.DataFrame, value: str, ylabel: str, title: str, subtitle: str = None
):
    figure, ax = plt.subplots(figsize=(10, 7.5))
    plotted = summary.replace([np.inf, -np.inf], np.nan).dropna(subset=["foundation", value, "lower", "upper"])
    below = np.clip(plotted[value] - plotted["lower"], 0, None)
    above = np.clip(plotted["upper"] - plotted[value], 0, None)
    ax.errorbar(
        plotted["foundation"],
        plotted[value],
        yerr=np.vstack([below, above]),
        fmt="o",
        color="#4D4D4D",
        ecolor="#BFBFBF",
        elinewidth=0.6,
        markersize=3.5,
        alpha=0.8,
    )
    ax.set_xlabel("Median foundation (BC negative, AD positive)")
    ax.set_ylabel(ylabel)
    ax.set_title(title if subtitle is None else "{}\n{}".format(title, subtitle))
    figure.tight_layout()
    return figure


def plot_by_size(summary: pd.DataFrame):
    figure, ax = plt.subplots(figsize=(10, 7.5))
    for name in SIZE_ORDER:
        part = summary[summary["size"] == name]
        ax.scatter(part["foundation"], part["mean"], s=16, alpha=0.7, color=SIZE_COLOURS[name], label=name)
        if len(part) >= 2:
            slope, intercept = np.polyfit(part["foundation"], part["mean"], 1)
            grid = np.linspace(part["foundation"].min(), part["foundation"].max(), 50)
            ax.plot(grid, intercept + slope * grid, color=SIZE_COLOURS[name], linewidth=2)
    ax.set_xlabel("Median foundation (BC negative, AD positive)")
    ax.set_ylabel("Mean persistence (years)")
    ax.set_title("Imputed persistence by size class")
    ax.legend(frameon=False, title="Size")
    figure.tight_layout()
    return figure


def plot_correlation(summary: pd.DataFrame):
    figure, ax = plt.subplots(figsize=(10, 7.5))
    for name in SIZE_ORDER:
        part = summary[summary["size"] == name]
        ax.scatter(part["foundation"], part["mean"], s=16, alpha=0.75, color=SIZE_COLOURS[name], label=name)
        if len(part) < 3:
            continue
        slope, intercept = np.polyfit(part["foundation"], part["mean"], 1)
        grid = np.linspace(part["foundation"].min(), part["foundation"].max(), 50)
        ax.plot(grid, intercept + slope * grid, color=SIZE_COLOURS[name], linewidth=2)
        coefficient, p_value = stats.pearsonr(part["foundation"], part["mean"])
        print("  {}: r={:.2f}, p={:.3g}, n={}".format(name, coefficient, p_value, len(part)))
    ax.set_xlabel("Median foundation (BC negative, AD positive)")
    ax.set_ylabel("Mean persistence (years)")
    ax.set_title("Persistence and foundation date")
    ax.legend(frameon=False)
    figure.tight_layout()
    return figure


def plot_persistence_counts(samples: pd.DataFrame, iteration: int):
    """One Monte Carlo iteration as counts, as in the R runner."""
    values = samples.loc[samples["iteration"] == iteration, "persistence"].to_numpy(dtype=float)
    values = values[np.isfinite(values)]
    figure, ax = plt.subplots(figsize=(10, 7.5))
    ax.hist(values, bins=ggplot_histogram_breaks(values), color="#CDAA7D", edgecolor="white", alpha=0.9)
    ax.set_xlabel("Persistence (years)")
    ax.set_ylabel("Count")
    ax.set_title("Site persistence, iteration {}".format(iteration))
    figure.tight_layout()
    return figure


def plot_parameter_draw(coefficients: pd.Series, seed: int = 1234):
    """1,000 draws from the mean generalised gamma, with no site data."""
    distribution = _prentice_gengamma(
        float(coefficients["mu"]), float(coefficients["sigma"]), float(coefficients["Q"])
    )
    figure, ax = plt.subplots(figsize=(10, 7.5))
    draws = distribution.rvs(1000, random_state=seed)
    ax.hist(
        draws,
        bins=r_histogram_breaks(draws),
        color="#D0D0D0",
        edgecolor="white",
    )
    ax.set_xlabel("Persistence (years)")
    ax.set_ylabel("Frequency")
    ax.set_title("Draws from the mean generalised gamma")
    figure.tight_layout()
    return figure


def plot_boxplot(summary: pd.DataFrame):
    figure, ax = plt.subplots(figsize=(10, 7.5))
    parts = [summary.loc[summary["size"] == name, "mean"].dropna() for name in SIZE_ORDER]
    drawn = ax.boxplot(parts, patch_artist=True)
    ax.set_xticks(range(1, len(SIZE_ORDER) + 1))
    ax.set_xticklabels(list(SIZE_ORDER))
    for box, name in zip(drawn["boxes"], SIZE_ORDER):
        box.set_facecolor(SIZE_COLOURS[name])
        box.set_alpha(0.6)
    ax.set_ylabel("Mean persistence (years)")
    ax.set_title("Imputed persistence by size class")
    figure.tight_layout()
    return figure


def run(settlements: mc.Settlements, destination: Path, sample_terminus: float, persistence_terminus: float, label) -> None:
    destination.mkdir(parents=True, exist_ok=True)
    print("Sampling {} settlements, {} iterations".format(settlements.n_settlements, ITERATIONS))
    trajectories = mc.sample_trajectories(
        settlements, ITERATIONS, terminus=sample_terminus, seed=SEED
    )
    samples = mc.mark_censoring(mc.persistence_samples(trajectories, terminus=persistence_terminus))
    print(
        "Persistence rows: {}; censored: {}; observed ends: {}".format(
            len(samples),
            int((samples["event"] == 0).sum()),
            int((samples["event"] == 1).sum()),
        )
    )

    print("Fitting a generalised gamma to each iteration")
    fitted = mc.fit_persistence(samples, distribution="gengamma")
    coefficients = fitted[["mu", "sigma", "Q"]].mean()
    print(
        "Mean coefficients: mu={:.4g} sigma={:.4g} Q={:.4g} ({} iterations)".format(
            coefficients["mu"], coefficients["sigma"], coefficients["Q"], len(fitted)
        )
    )
    fitted.to_csv(destination / "gengamma-coefficients.csv", index=False)

    save(plot_parameter_draw(coefficients), destination, "gengamma-parameter-draw")
    save(plot_persistence_counts(samples, iteration=123), destination, "persistence-iteration-123")
    histogram = mc.plot_longevity(samples)
    histogram.axes[0].set_title("Site persistence across draws")
    save(histogram, destination, "persistence-histogram")
    save(
        mc.plot_longevity(samples, coefficients=coefficients),
        destination,
        "persistence-gengamma-overlay",
    )

    # longevity.R replaces every draw of a site that was censored in any
    # iteration, and leaves norm_pers at the value from before imputation.
    imputed = samples.copy()
    still_alive = imputed.loc[imputed["event"] == 0, "id"].unique()
    replace = imputed["id"].isin(still_alive)
    imputed.loc[replace, "persistence"] = mc.impute_lifespan(
        imputed.loc[replace, "persistence"].to_numpy(),
        coefficients["mu"],
        coefficients["sigma"],
        coefficients["Q"],
        seed=1234,
    )
    summary = mc.summarize_persistence(imputed)
    summary["foundation"] = mc.bcad_to_bp(summary["start"], bc_to_bp=False)
    if label is not None:
        summary["size"] = pd.Categorical(summary["id"].map(label), categories=SIZE_ORDER, ordered=True)
    summary.to_csv(destination / "persistence-by-site.csv", index=False)

    if label is not None:
        by_size = (
            summary.groupby("size", observed=True)
            .agg(
                mean_persistence=("mean", "mean"),
                median_persistence=("mean", "median"),
                q975=("mean", lambda values: values.quantile(0.975)),
                q025=("mean", lambda values: values.quantile(0.025)),
                n=("mean", "size"),
            )
            .reset_index()
        )
        by_size.to_csv(destination / "persistence-by-size.csv", index=False)
        print(by_size.to_string(index=False))

    save(
        plot_against_foundation(
            summary.assign(lower=summary["q025"], upper=summary["q975"]),
            "mean",
            "Mean persistence (years)",
            "Imputed persistence against foundation date",
        ),
        destination,
        "persistence-vs-foundation",
    )
    save(
        plot_against_foundation(
            summary.assign(lower=summary["q025_n"], upper=summary["q975_n"]),
            "mean_n",
            "Mean share of available time",
            "Normalised persistence against foundation date",
            subtitle="Normalised values are from before lifespan imputation",
        ),
        destination,
        "norm-persistence-vs-foundation",
    )
    if label is not None:
        save(plot_by_size(summary), destination, "persistence-vs-foundation-by-size")
        save(plot_correlation(summary), destination, "persistence-correlation-by-size")
        save(plot_boxplot(summary), destination, "persistence-boxplot-by-size")
    print("Longevity figures written to {}".format(destination.relative_to(ROOT)))


def main() -> None:
    parser = argparse.ArgumentParser(description="Longevity figures for an example dataset.")
    parser.add_argument("--dataset", choices=("yautepec", "viabundus"), default="yautepec")
    args = parser.parse_args()
    if args.dataset == "yautepec":
        source = DATA / "yautepec.csv"
        run(read_yautepec(source), FIGURES / "yautepec" / "longevity-python", 510, 431, size_labels(source))
    else:
        settlements = census_phases(read_viabundus(DATA / "viabundus.csv"))
        run(settlements, FIGURES / "viabundus" / "longevity-python", 1625, 1675, None)


if __name__ == "__main__":
    main()
