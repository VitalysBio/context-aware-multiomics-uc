from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from multiomics_uc.paths import get_path_from_root


def load_expression() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_expression_log2_cpm_filtered.csv"
    )
    return pd.read_csv(path, index_col=0)


def load_clusters() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)[["patient_id", "rnaseq_cluster"]]


def align_expression_and_clusters(
    expression: pd.DataFrame,
    clusters: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    clusters = clusters.set_index("patient_id")

    common_patients = expression.index.intersection(clusters.index)

    expression = expression.loc[common_patients]
    clusters = clusters.loc[common_patients]

    return expression, clusters


def compute_cluster_markers(
    expression: pd.DataFrame,
    clusters: pd.DataFrame,
    cluster_label: str,
) -> pd.DataFrame:
    in_cluster = clusters["rnaseq_cluster"].astype(str) == str(cluster_label)

    expr_cluster = expression.loc[in_cluster]
    expr_rest = expression.loc[~in_cluster]

    mean_cluster = expr_cluster.mean(axis=0)
    mean_rest = expr_rest.mean(axis=0)

    log2_fc = mean_cluster - mean_rest

    statistic, p_values = ttest_ind(
        expr_cluster,
        expr_rest,
        axis=0,
        equal_var=False,
        nan_policy="omit",
    )

    _, adj_p_values, _, _ = multipletests(
        p_values,
        alpha=0.05,
        method="fdr_bh",
    )

    result = pd.DataFrame(
        {
            "gene_id": expression.columns,
            "mean_cluster": mean_cluster.values,
            "mean_rest": mean_rest.values,
            "log2_fc": log2_fc.values,
            "t_statistic": statistic,
            "p_value": p_values,
            "adj_p_value": adj_p_values,
        }
    )

    result["abs_log2_fc"] = result["log2_fc"].abs()

    result = result.sort_values(
        by=["adj_p_value", "abs_log2_fc"],
        ascending=[True, False],
    )

    return result


def save_cluster_markers(markers: pd.DataFrame, cluster_label: str) -> None:
    out_dir = get_path_from_root("reports/tables/cluster_markers")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"rnaseq_cluster_{cluster_label}_markers_vs_rest.csv"
    markers.to_csv(out_file, index=False)

    print(f"[OK] Saved markers for cluster {cluster_label}: {out_file}")


def summarize(markers: pd.DataFrame, cluster_label: str) -> None:
    significant = markers[
        (markers["adj_p_value"] < 0.05) &
        (markers["abs_log2_fc"] >= 1.0)
    ]

    up = significant[significant["log2_fc"] > 0]
    down = significant[significant["log2_fc"] < 0]

    print(f"\n=== Cluster {cluster_label} markers ===")
    print(f"Significant markers, FDR < 0.05 and |log2FC| >= 1: {len(significant)}")
    print(f"Up in cluster: {len(up)}")
    print(f"Down in cluster: {len(down)}")

    print("\nTop upregulated genes:")
    print(
        up.sort_values(["adj_p_value", "log2_fc"], ascending=[True, False])
        .head(10)[["gene_id", "log2_fc", "adj_p_value"]]
        .to_string(index=False)
    )


def main() -> int:
    try:
        expression = load_expression()
        clusters = load_clusters()

        expression, clusters = align_expression_and_clusters(
            expression,
            clusters,
        )

        cluster_labels = sorted(clusters["rnaseq_cluster"].astype(str).unique())

        print("\n Differential Expression by Cluster")
        print(f"Patients: {expression.shape[0]}")
        print(f"Genes tested: {expression.shape[1]}")
        print(f"Clusters: {cluster_labels}")

        for cluster_label in cluster_labels:
            markers = compute_cluster_markers(
                expression,
                clusters,
                cluster_label,
            )
            summarize(markers, cluster_label)
            save_cluster_markers(markers, cluster_label)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())