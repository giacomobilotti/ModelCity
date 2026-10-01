"""Chronological Monte Carlo sampling for every example dataset.

Runs ``modelcity.sample_trajectories`` with the settings of the R sampling
scripts and writes per-settlement growth and trajectory figures, period growth,
persistence, and the matching summary tables in the column layout of the R
output.  Yautepec mirrors ``scripts/r/run-sampling-plots.R``; Viabundus samples
every city and plots a selection.

Viabundus records population at census years rather than over phases, so each
census is first read as the period ``as_settlements`` infers around it
(1300 covers 1250-1350, and so on) and zero-population rows are dropped, since
the sampling treats every row as an occupied phase.

    python scripts/run_sampling.py
    python scripts/run_sampling.py --dataset viabundus
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

import modelcity as mc  # noqa: E402
from run_unified import DATA, FIGURES, ROOT, read_viabundus, read_yautepec  # noqa: E402

ITERATIONS = 500
SEED = 42


def census_phases(settlements: mc.Settlements) -> mc.Settlements:
    """Snapshot records as interval phases spanning their inferred periods.

    Rows with no population are dropped: a phase in the sampling is a period
    the settlement occupied, so a later first census reads as a later
    foundation and an earlier last census as an earlier abandonment.  The
    sampling reads each settlement's phases in row order, so rows are sorted
    oldest first.
    """
    data = settlements.data
    occupied = (
        data[data["size"] > 0]
        .assign(_order=lambda frame: pd.factorize(frame["id"])[0])
        .sort_values(["_order", "t_start"], kind="stable")
        .drop(columns="_order")
        .reset_index(drop=True)
    )
    time = mc.TimeSpec(time_scale=settlements.time.time_scale, record_type="interval")
    return mc.Settlements(occupied, settlements.size, time, settlements.roles)


def largest(settlements: mc.Settlements, ids, n: int) -> List[str]:
    sizes = settlements.data[settlements.data["id"].isin(ids)].groupby("id")["size"].max()
    return sizes.sort_values(ascending=False).index[:n].tolist()


def yautepec_focus(settlements: mc.Settlements) -> List[str]:
    return ["160", "200", "126", "133", "148", "334", "341"]


def viabundus_focus(settlements: mc.Settlements) -> List[str]:
    """The largest cities seen at every census, plus the largest that appear late or vanish early."""
    data = settlements.data
    first = data.groupby("id")["t_start"].min()
    last = data.groupby("id")["t_end"].max()
    full = first.index[(first == first.min()) & (last == last.max())]
    late = first.index[first > first.min()]
    early_end = last.index[last < last.max()]
    return largest(settlements, full, 5) + largest(settlements, late, 2) + largest(settlements, early_end, 1)


@dataclass(frozen=True)
class Run:
    load: Callable[[], mc.Settlements]
    folder: str
    focus: Callable[[mc.Settlements], List[str]]
    #: Settlements whose sampled end falls before this year end at size zero.
    terminus: float
    #: End of the record for persistence.
    persistence_terminus: float
    #: Sample only the focus settlements, as the R runner does, rather than all.
    focus_only: bool = False
    #: Settlements resampled with exponential dates and normal peaks.
    expnorm: Tuple[str, ...] = ()


#: Terminus values follow the R examples: the start of the last phase for
#: sampling and its end for persistence, in each dataset's own time scale.
#: Yautepec repeats ``scripts/r/run-sampling-plots.R`` so the two folders
#: can be compared file by file.
RUNS: Dict[str, Run] = {
    "yautepec": Run(
        load=lambda: read_yautepec(DATA / "yautepec.csv"),
        folder="yautepec/sampling-python",
        focus=yautepec_focus,
        terminus=510,
        persistence_terminus=431,
        focus_only=True,
        expnorm=("160", "200"),
    ),
    "viabundus": Run(
        load=lambda: census_phases(read_viabundus(DATA / "viabundus.csv")),
        folder="viabundus/sampling-python",
        focus=viabundus_focus,
        terminus=1625,
        persistence_terminus=1675,
    ),
}


def run_dataset(key: str) -> None:
    run = RUNS[key]
    settlements = run.load()
    destination = FIGURES / run.folder
    focus = run.focus(settlements)
    names = settlements.data.drop_duplicates("id").set_index("id")["name"]

    print("\n{}".format(key))
    print(settlements)
    print("  focus: {}".format(", ".join("{} ({})".format(i, names[i]) for i in focus)))

    sampled = only(settlements, focus) if run.focus_only else settlements
    trajectories = mc.sample_trajectories(sampled, ITERATIONS, terminus=run.terminus, seed=SEED)
    focused = trajectories[trajectories["id"].isin(focus)]

    def save(figure, name: str) -> None:
        path = mc.save_figure(figure, destination / "{}.png".format(name))
        plt.close(figure)
        print("  wrote {}".format(Path(path).relative_to(ROOT)))

    # As in R, single-phase settlements get no plots: their growth has one point.
    phases = settlements.data.groupby("id").size()
    for site in focus:
        if phases[site] < 2:
            print("  skipped {} (single phase)".format(site))
            continue
        save(mc.plot_site_growth(focused, site), "cagr-site-{}".format(site))
        save(mc.plot_trajectories(focused, site, value="growth_cagr"), "spaghetti-cagr-site-{}".format(site))
        save(
            mc.plot_trajectories(focused, site, value="size"),
            "spaghetti-{}-site-{}".format(settlements.size.size_type, site),
        )
    save(mc.plot_period_growth(trajectories), "cagr-by-period")
    save(mc.plot_persistence(focused, terminus=run.persistence_terminus), "persistence-by-site")

    if run.expnorm:
        variant = mc.sample_trajectories(
            only(settlements, run.expnorm),
            ITERATIONS,
            start_method="exp",
            end_method="exp",
            peak_method="norm",
            perc=5,
            peak_sd_perc=0.1,
            terminus=run.terminus,
            seed=SEED,
        )
        for site in run.expnorm:
            if phases[site] < 2:
                print("  skipped exp/norm {} (single phase)".format(site))
                continue
            save(mc.plot_site_growth(variant, site), "cagr-expnorm-site-{}".format(site))
            save(
                mc.plot_trajectories(variant, site, value="growth_cagr"),
                "spaghetti-expnorm-cagr-site-{}".format(site),
            )

    tables = {
        "summary-by-site": mc.site_summary(trajectories),
        "summary-by-period": mc.period_summary(trajectories),
        "persistence-by-site": mc.persistence_summary(trajectories, terminus=run.persistence_terminus),
    }
    for name, table in tables.items():
        path = destination / "{}.csv".format(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        r_layout(table, names, settlements.size.size_type).to_csv(path, index=False)
        print("  wrote {}".format(path.relative_to(ROOT)))


def only(settlements: mc.Settlements, ids) -> mc.Settlements:
    data = settlements.data
    return settlements.with_data(data[data["id"].isin(ids)].reset_index(drop=True))


def r_layout(table: pd.DataFrame, names: pd.Series, size_type: str) -> pd.DataFrame:
    """A summary table with the columns of the R ``summary_*`` and ``persistence`` output.

    Phase bounds become ``start`` and ``end`` in years BP and the size column
    takes the size type's name (R calls it ``area``).  Settlement names are
    kept after ``id`` when they say more than the identifier.
    """
    frame = table.drop(columns=["label", "year_median"], errors="ignore")
    if "t_start" in frame.columns:
        frame["t_start"] = mc.to_bp(frame["t_start"])
        frame["t_end"] = mc.to_bp(frame["t_end"])
    frame = frame.rename(columns={"t_start": "start", "t_end": "end", "size": size_type})
    if "id" in frame.columns:
        labels = frame["id"].map(names)
        if (labels.astype(str) != frame["id"].astype(str)).any():
            frame.insert(1, "name", labels)
    return frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        choices=sorted(RUNS),
        action="append",
        help="Dataset to sample; repeat the flag for several, default is all.",
    )
    arguments = parser.parse_args()

    for key in arguments.dataset or sorted(RUNS):
        run_dataset(key)

    print("\nDone.")


if __name__ == "__main__":
    main()
