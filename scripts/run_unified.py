"""Regenerate the core figures for every dataset through the shared API.

The point of this script is that the two datasets differ only in how they are
loaded.  Once each is a ``Settlements`` object, the same seven calls produce the
same seven figures, with areas in hectares for one and headcounts for the other.

    python3 scripts/run_unified.py
    python3 scripts/run_unified.py --dataset yautepec
    python3 scripts/run_unified.py --years BCAD
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, List, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402

import modelcity as mc  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "examples"
FIGURES = ROOT / "figures"

#: Each dataset is a loader plus the directory its figures belong in.
DATASETS: Dict[str, Tuple[Callable[[], "mc.Settlements"], str]] = {
    "yautepec": (lambda: mc.read_yautepec(DATA / "yautepec.csv"), "yautepec/unified"),
    "viabundus": (lambda: mc.read_viabundus(DATA / "viabundus.csv"), "viabundus/unified"),
}

#: Figure name to plotting function, identical for every dataset.
FIGURE_SET: List[Tuple[str, Callable]] = [
    ("01-duration", mc.plot_duration),
    ("02-foundations-abandonments", mc.plot_events),
    ("03-churn", mc.plot_churn),
    ("04-size-totals", mc.plot_size_totals),
    ("05-size-summary", mc.plot_size_summary),
    ("06-rank-size", mc.plot_rank_size),
    ("07-urban-share", mc.plot_urban_share),
]

#: Metric tables written alongside the figures.
TABLE_SET: List[Tuple[str, Callable]] = [
    ("tbl-duration", mc.settlement_duration),
    ("tbl-size-totals", mc.size_totals),
    ("tbl-size-summary", mc.size_summary),
    ("tbl-urban-share", mc.urban_share),
]


def run_dataset(key: str) -> None:
    loader, folder = DATASETS[key]
    settlements = loader()
    destination = FIGURES / folder

    print("\n{}".format(key))
    print(settlements)

    for name, plot in FIGURE_SET:
        figure = plot(settlements)
        path = mc.save_figure(figure, destination / "{}.png".format(name))
        plt.close(figure)
        print("  wrote {}".format(Path(path).relative_to(ROOT)))

    for name, metric in TABLE_SET:
        table = metric(settlements)
        path = destination / "{}.csv".format(name)
        path.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(path, index=False)
        print("  wrote {}".format(path.relative_to(ROOT)))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        choices=sorted(DATASETS),
        action="append",
        help="Dataset to render; repeat the flag for several, default is all.",
    )
    parser.add_argument(
        "--years",
        choices=mc.YEAR_DISPLAYS,
        default=mc.get_year_display(),
        help="Scale for axis ticks and period labels; affects labels only.",
    )
    arguments = parser.parse_args()

    mc.set_year_display(arguments.years)
    for key in arguments.dataset or sorted(DATASETS):
        run_dataset(key)

    print("\nDone.")


if __name__ == "__main__":
    main()
