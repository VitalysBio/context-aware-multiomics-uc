from __future__ import annotations

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from multiomics_uc.paths import get_path_from_root


PROGRAMS = {
    "immune_infiltration_percentile": "Immune\ninfiltration",
    "angiogenic_stroma_percentile": "Angiogenic\nstroma",
    "proliferation_percentile": "Proliferation",
    "motility_cytoskeleton_percentile": "Motility /\ncytoskeleton",
    "immune_activation_percentile": "IFN immune\nactivation",
    "invasion_signaling_percentile": "Invasion\nsignaling",
    "differentiation_percentile": "Differentiation",
}


def load_percentiles() -> pd.DataFrame:
    path = get_path_from_root("data/processed/program_percentiles.csv")
    return pd.read_csv(path)


def get_patient_row(df: pd.DataFrame, patient_id: str) -> pd.Series:
    subset = df[df["patient_id"] == patient_id]

    if subset.empty:
        raise ValueError(f"Patient not found: {patient_id}")

    return subset.iloc[0]


def plot_radar(patient_id: str, row: pd.Series) -> None:
    labels = list(PROGRAMS.values())
    cols = list(PROGRAMS.keys())

    values = [float(row[col]) for col in cols]

    values += values[:1]

    angles = np.linspace(
        0,
        2 * np.pi,
        len(labels),
        endpoint=False,
    ).tolist()

    angles += angles[:1]

    fig, ax = plt.subplots(
        figsize=(8, 8),
        subplot_kw={"polar": True},
    )

    ax.plot(
        angles,
        values,
        linewidth=2,
    )

    ax.fill(
        angles,
        values,
        alpha=0.25,
    )

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(labels, fontsize=10)

    ax.set_ylim(0, 100)

    ax.set_yticks([20, 40, 60, 80, 100])
    ax.set_yticklabels(["20", "40", "60", "80", "100"])

    ax.set_title(
        f"Biological Program Profile\n{patient_id}",
        fontsize=14,
        weight="bold",
        pad=25,
    )

    figures_dir = get_path_from_root("reports/figures/patient_reports")
    figures_dir.mkdir(parents=True, exist_ok=True)

    out_file = figures_dir / f"{patient_id}_biological_program_radar.png"

    fig.savefig(
        out_file,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"[OK] Saved radar plot: {out_file}")


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--patient-id",
        required=True,
        help="TCGA patient ID, e.g. TCGA-2F-A9KO",
    )

    args = parser.parse_args()

    try:
        percentiles = load_percentiles()
        row = get_patient_row(percentiles, args.patient_id)

        plot_radar(
            args.patient_id,
            row,
        )

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())