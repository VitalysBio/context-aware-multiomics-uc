from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from multiomics_uc.paths import get_path_from_root


PROGRAMS = [
    {
        "pc": "PC1",
        "title": "Immune-stromal\nmicroenvironment",
        "biology": "Macrophages\nComplement\nECM / EMT",
        "cluster": "Cluster 3",
        "role": "RSF important",
    },
    {
        "pc": "PC3",
        "title": "Proliferation\nprogram",
        "biology": "Cell cycle\nMitosis\nChromosome segregation",
        "cluster": "Cluster 2",
        "role": "RSF important",
    },
    {
        "pc": "PC6",
        "title": "Immune activation\nvs EMT",
        "biology": "IFNγ / IFNα\nCytokines\nECM / fibrosis",
        "cluster": "Clusters 3 / 0",
        "role": "RSF important",
    },
    {
        "pc": "PC8",
        "title": "Invasion and\nRTK signaling",
        "biology": "RHO / RAC\nMET motility\nRTK signaling",
        "cluster": "Cluster 3",
        "role": "Cox risk",
    },
    {
        "pc": "PC9",
        "title": "Differentiation\nprogram",
        "biology": "Cilium assembly\nEpithelial organization\nCell projections",
        "cluster": "Cluster 1",
        "role": "Cox protective",
    },
]


def add_box(ax, x, y, width, height, text, fontsize=9, weight="normal"):
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.04,rounding_size=0.08",
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
        weight=weight,
    )


def add_arrow(ax, start, end):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="->",
        mutation_scale=14,
        linewidth=1.1,
        color="black",
    )
    ax.add_patch(arrow)


def main() -> int:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(15, 8))
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 8)
    ax.axis("off")

    ax.text(
        7.5,
        7.55,
        "Latent Biological Programs Driving TCGA-BLCA Heterogeneity",
        ha="center",
        va="center",
        fontsize=16,
        weight="bold",
    )

    ax.text(
        7.5,
        7.15,
        "RNA-seq PCA components interpreted through pathway enrichment, cluster association, and survival modeling",
        ha="center",
        va="center",
        fontsize=10,
    )

    x_positions = [0.6, 3.45, 6.3, 9.15, 12.0]

    for x, program in zip(x_positions, PROGRAMS):
        add_box(
            ax,
            x,
            5.75,
            2.25,
            0.65,
            program["pc"],
            fontsize=13,
            weight="bold",
        )

        add_box(
            ax,
            x,
            4.75,
            2.25,
            0.8,
            program["title"],
            fontsize=10,
            weight="bold",
        )

        add_box(
            ax,
            x,
            3.35,
            2.25,
            1.15,
            program["biology"],
            fontsize=9,
        )

        add_box(
            ax,
            x,
            2.25,
            2.25,
            0.75,
            f"Associated with\n{program['cluster']}",
            fontsize=8,
        )

        add_box(
            ax,
            x,
            1.25,
            2.25,
            0.65,
            program["role"],
            fontsize=8,
        )

        add_arrow(ax, (x + 1.125, 5.75), (x + 1.125, 5.55))
        add_arrow(ax, (x + 1.125, 4.75), (x + 1.125, 4.50))
        add_arrow(ax, (x + 1.125, 3.35), (x + 1.125, 3.00))
        add_arrow(ax, (x + 1.125, 2.25), (x + 1.125, 1.90))

    add_box(
        ax,
        2.2,
        0.25,
        10.6,
        0.65,
        "Interpretation: TCGA-BLCA prognosis is shaped by continuous latent programs \n of immune-stromal infiltration, proliferation, immune activation, invasion, and differentiation.",
        fontsize=9,
    )

    out_file = figures_dir / "latent_programs_summary.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Saved latent programs summary figure: {out_file}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())