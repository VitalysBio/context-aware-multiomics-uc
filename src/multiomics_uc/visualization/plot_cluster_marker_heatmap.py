from __future__ import annotations

import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import zscore

from multiomics_uc.paths import get_path_from_root


TOP_GENES_PER_CLUSTER = 15


def load_expression() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_expression_log2_cpm_filtered.csv"
    )
    return pd.read_csv(path, index_col=0)


def load_clusters() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_rnaseq_clusters.csv"
    )
    return pd.read_csv(path)[["patient_id", "rnaseq_cluster"]]


def load_marker_file(cluster_label: str) -> pd.DataFrame:
    path = get_path_from_root(
        f"reports/tables/cluster_markers_annotated/"
        f"rnaseq_cluster_{cluster_label}_markers_vs_rest_annotated.csv"
    )
    return pd.read_csv(path)


def get_top_marker_genes() -> list[str]:
    genes = []

    for cluster in ["0", "1", "2", "3"]:
        markers = load_marker_file(cluster)

        markers = markers[
            (markers["adj_p_value"] < 0.05)
            & (markers["log2_fc"] > 1)
            & (markers["gene_name"].notna())
            & (markers["gene_type"] == "protein_coding")
        ].copy()

        top = (
            markers
            .sort_values(["adj_p_value", "log2_fc"], ascending=[True, False])
            .head(TOP_GENES_PER_CLUSTER)
        )

        genes.extend(top["gene_name"].tolist())

    genes = list(dict.fromkeys(genes))

    return genes


def build_gene_symbol_mapping() -> dict:
    annotation = pd.read_csv(
        get_path_from_root("data/processed/tcga_blca_gene_annotation_table.csv")
    )

    mapping = (
        annotation[["gene_id", "gene_name"]]
        .drop_duplicates()
        .set_index("gene_name")["gene_id"]
        .to_dict()
    )

    return mapping


def build_heatmap_matrix(
    expression: pd.DataFrame,
    clusters: pd.DataFrame,
    selected_genes: list[str],
    mapping: dict,
) -> tuple[pd.DataFrame, pd.Series]:
    gene_ids = [
        mapping[g]
        for g in selected_genes
        if g in mapping
    ]

    gene_labels = [
        g
        for g in selected_genes
        if g in mapping
    ]

    matrix = expression[gene_ids].copy()

    matrix.columns = gene_labels

    clusters_sorted = clusters.sort_values("rnaseq_cluster")
    
    matrix = matrix.loc[clusters_sorted["patient_id"]]
    
    cluster_labels = (
    clusters_sorted
    .set_index("patient_id")
    .loc[matrix.index]["rnaseq_cluster"]
    )

    matrix_z = matrix.apply(zscore, axis=0)

    return matrix_z.T, cluster_labels


def plot_heatmap(
    matrix_z: pd.DataFrame,
    cluster_labels: pd.Series,
) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    lut = {
        "0": "#8dd3c7",
        "1": "#ffffb3",
        "2": "#fb8072",
        "3": "#80b1d3",
    }

    col_colors = cluster_labels.astype(str).map(lut)

    g = sns.clustermap(
        matrix_z,
        col_cluster=False,
        row_cluster=True,
        col_colors=col_colors,
        cmap="vlag",
        figsize=(16, 12),
        xticklabels=False,
        yticklabels=True,
        center=0,
    )

    g.fig.suptitle(
        "Top transcriptomic marker genes across RNA-seq clusters",
        y=1.02,
        fontsize=16,
    )

    out_file = figures_dir / "tcga_blca_cluster_marker_heatmap.png"

    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Saved heatmap: {out_file}")


def summarize(matrix_z: pd.DataFrame) -> None:
    print("\n Cluster Marker Heatmap Summary")
    print(f"Genes plotted: {matrix_z.shape[0]}")
    print(f"Patients plotted: {matrix_z.shape[1]}")


def main() -> int:
    try:
        expression = load_expression()
        clusters = load_clusters()

        selected_genes = get_top_marker_genes()
        mapping = build_gene_symbol_mapping()

        matrix_z, cluster_labels = build_heatmap_matrix(
            expression,
            clusters,
            selected_genes,
            mapping,
        )

        summarize(matrix_z)
        plot_heatmap(matrix_z, cluster_labels)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())