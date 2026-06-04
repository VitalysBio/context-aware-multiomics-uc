from __future__ import annotations

import sys

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from scipy.stats import kruskal

from multiomics_uc.paths import get_path_from_root


PCS = ["PC1", "PC2", "PC3", "PC4", "PC6", "PC8", "PC9"]


def load_embeddings() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_rnaseq_pca_umap_embeddings.csv"
    )
    return pd.read_csv(path)


def load_clusters() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_rnaseq_clusters.csv"
    )
    return pd.read_csv(path)[["patient_id", "rnaseq_cluster"]]


def merge_embeddings_clusters(
    embeddings: pd.DataFrame,
    clusters: pd.DataFrame,
) -> pd.DataFrame:
    df = embeddings.merge(
        clusters,
        on="patient_id",
        how="inner",
        validate="one_to_one",
    )

    df["rnaseq_cluster"] = df["rnaseq_cluster"].astype(str)

    return df


def kruskal_test(df: pd.DataFrame, pc: str) -> dict:
    groups = []

    for cluster in sorted(df["rnaseq_cluster"].unique()):
        values = df.loc[df["rnaseq_cluster"] == cluster, pc].dropna()
        groups.append(values)

    statistic, pvalue = kruskal(*groups)

    return {
        "pc": pc,
        "kruskal_statistic": statistic,
        "p_value": pvalue,
    }


def plot_boxplot(df: pd.DataFrame, pc: str) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 5))

    sns.boxplot(
        data=df,
        x="rnaseq_cluster",
        y=pc,
        order=sorted(df["rnaseq_cluster"].unique()),
    )

    sns.stripplot(
        data=df,
        x="rnaseq_cluster",
        y=pc,
        order=sorted(df["rnaseq_cluster"].unique()),
        alpha=0.3,
        size=2,
    )

    plt.title(f"{pc} by RNA-seq Cluster")
    plt.xlabel("RNA-seq Cluster")
    plt.ylabel(pc)

    plt.tight_layout()

    out_file = figures_dir / f"{pc}_by_cluster_boxplot.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"[OK] Saved: {out_file}")


def summarize(df: pd.DataFrame, pc: str) -> None:
    print(f"\n=== {pc} by Cluster ===")

    summary = (
        df.groupby("rnaseq_cluster")[pc]
        .agg(["count", "mean", "median", "std"])
        .round(3)
    )

    print(summary)


def save_stats(results: list[dict]) -> None:
    out_dir = get_path_from_root(
        "reports/tables/pca_interpretation"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "pc_cluster_association_statistics.csv"

    pd.DataFrame(results).to_csv(
        out_file,
        index=False,
    )

    print(f"\n[OK] Saved statistics: {out_file}")


def main() -> int:
    try:
        embeddings = load_embeddings()
        clusters = load_clusters()

        df = merge_embeddings_clusters(
            embeddings,
            clusters,
        )

        print("\n PC Cluster Association Analysis")
        print(f"Patients: {df['patient_id'].nunique()}")
        print(f"Clusters: {sorted(df['rnaseq_cluster'].unique())}")

        results = []

        for pc in PCS:
            if pc not in df.columns:
                raise ValueError(f"{pc} not found in embeddings table.")

            summarize(df, pc)

            stats = kruskal_test(df, pc)
            results.append(stats)

            print(
                f"\nKruskal-Wallis statistic: "
                f"{stats['kruskal_statistic']:.4f}"
            )

            print(
                f"Kruskal-Wallis p-value: "
                f"{stats['p_value']:.4e}"
            )

            plot_boxplot(df, pc)

        save_stats(results)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())