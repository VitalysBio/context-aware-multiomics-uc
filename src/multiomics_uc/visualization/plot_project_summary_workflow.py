from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from multiomics_uc.paths import get_path_from_root


def add_box(ax, xy, text, width=2.4, height=0.75, fontsize=9):
    x, y = xy
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.03,rounding_size=0.08",
        linewidth=1.2,
        edgecolor="black",
        facecolor="white",
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        wrap=True,
    )


def add_arrow(ax, start, end):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="->",
        mutation_scale=14,
        linewidth=1.2,
        color="black",
    )
    ax.add_patch(arrow)


def main() -> int:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(14, 9))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 9)
    ax.axis("off")

    add_box(ax, (0.7, 7.3), "TCGA-BLCA\nclinical + survival", 2.4)
    add_box(ax, (3.9, 7.3), "RNA-seq\n406 patients", 2.4)
    add_box(ax, (7.1, 7.3), "Somatic mutations\n406 patients", 2.4)
    add_box(ax, (10.3, 7.3), "Gene-level CNV\n197 patients", 2.4)

    add_box(ax, (3.9, 5.8), "PCA + UMAP\nRNA embeddings", 2.4)
    add_box(ax, (3.9, 4.4), "Transcriptomic\nclusters", 2.4)

    add_box(ax, (0.7, 3.0), "Cluster biology\nECM, luminal,\nbasal, immune", 2.6)
    add_box(ax, (3.9, 3.0), "Survival association\nKM + Cox", 2.4)
    add_box(ax, (7.1, 3.0), "Driver landscape\nTP53, RB1,\nKDM6A", 2.4)
    add_box(ax, (10.3, 3.0), "CNV architecture\nEGFR, MYC,\nCDKN2A", 2.4)

    add_box(ax, (2.4, 1.3), "Integrated patient\nfeature table", 2.8)
    add_box(ax, (5.8, 1.3), "Predictive benchmark\nCox + RSF", 2.8)
    add_box(ax, (9.2, 1.3), "Interpretable latent\nprograms\nPC1, PC3, PC6,\nPC8, PC9", 3.0)

    add_arrow(ax, (1.9, 7.3), (3.9, 6.2))
    add_arrow(ax, (5.1, 7.3), (5.1, 6.55))
    add_arrow(ax, (5.1, 5.8), (5.1, 5.15))
    add_arrow(ax, (5.1, 4.4), (5.1, 3.75))

    add_arrow(ax, (8.3, 7.3), (8.3, 3.75))
    add_arrow(ax, (11.5, 7.3), (11.5, 3.75))

    add_arrow(ax, (2.0, 3.0), (3.5, 2.05))
    add_arrow(ax, (5.1, 3.0), (4.0, 2.05))
    add_arrow(ax, (8.3, 3.0), (4.9, 2.05))
    add_arrow(ax, (11.5, 3.0), (4.9, 2.05))

    add_arrow(ax, (5.2, 1.68), (5.8, 1.68))
    add_arrow(ax, (8.6, 1.68), (9.2, 1.68))

    ax.text(
        7,
        8.55,
        "Context-aware Multi-Omics Patient Stratification and Survival Prediction",
        ha="center",
        fontsize=15,
        weight="bold",
    )

    ax.text(
        7,
        0.35,
        "Key result: biologically coherent transcriptomic tumor states, complementary mutation/CNV landscapes,\n"
        "and cross-validated clinical + RNA + mutation survival prediction.",
        ha="center",
        fontsize=10,
    )

    out_file = figures_dir / "project_summary_workflow.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Saved project workflow figure: {out_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())