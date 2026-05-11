from __future__ import annotations

import sys
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from multiomics_uc.paths import get_path_from_root


CLUSTER_LABELS = {
    "0": "Stromal / ECM-rich",
    "1": "Luminal-like metabolic",
    "2": "Basal / squamous inflammatory",
    "3": "Immune-inflamed IFNγ-high",
}


def clean_term(term: str) -> str:
    term = str(term)
    term = term.split(" R-HSA")[0]
    term = term.replace("(GO:", "\n(GO:")
    return "\n".join(textwrap.wrap(term, width=38))


def load_enrichment(cluster: str) -> pd.DataFrame:
    path = get_path_from_root(
        "reports/tables/pathway_enrichment",
        f"rnaseq_cluster_{cluster}_pathway_enrichment.csv",
    )
    return pd.read_csv(path)


def select_top_pathways(df: pd.DataFrame, n: int = 8) -> pd.DataFrame:
    df = df.copy()
    df = df[df["Adjusted P-value"].notna()]
    df["neg_log10_fdr"] = -np.log10(df["Adjusted P-value"])

    df = df.sort_values("Adjusted P-value").head(n)
    df["clean_term"] = df["Term"].apply(clean_term)

    return df.sort_values("neg_log10_fdr", ascending=True)


def plot_summary() -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    axes = axes.flatten()

    for ax, cluster in zip(axes, ["0", "1", "2", "3"]):
        enrichment = load_enrichment(cluster)
        top = select_top_pathways(enrichment, n=8)

        ax.barh(top["clean_term"], top["neg_log10_fdr"])
        ax.set_title(f"Cluster {cluster}: {CLUSTER_LABELS[cluster]}")
        ax.set_xlabel("-log10 adjusted p-value")
        ax.tick_params(axis="y", labelsize=8)

    fig.suptitle(
        "Top enriched pathways across RNA-seq transcriptomic clusters",
        fontsize=16,
        y=1.02,
    )

    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_pathway_enrichment_summary.png"
    fig.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[OK] Saved pathway summary figure: {out_file}")


def main() -> int:
    try:
        plot_summary()
        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())