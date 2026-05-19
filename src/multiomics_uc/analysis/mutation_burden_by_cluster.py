from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import kruskal

from multiomics_uc.paths import get_path_from_root


def load_clusters() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)[["patient_id", "rnaseq_cluster"]]


def load_burden() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_mutation_burden.csv")
    return pd.read_csv(path)


def merge_data(clusters: pd.DataFrame, burden: pd.DataFrame) -> pd.DataFrame:
    df = clusters.merge(burden, on="patient_id", how="left")
    df["mutation_burden"] = df["mutation_burden"].fillna(0)
    df["rnaseq_cluster"] = df["rnaseq_cluster"].astype(str)
    return df


def run_test(df: pd.DataFrame) -> tuple[float, float]:
    groups = [
        subset["mutation_burden"].dropna()
        for _, subset in df.groupby("rnaseq_cluster")
    ]
    statistic, p_value = kruskal(*groups)
    return statistic, p_value


def summarize(df: pd.DataFrame, statistic: float, p_value: float) -> None:
    print("\n Mutation Burden by RNA-seq Cluster")
    print(f"Patients: {df['patient_id'].nunique()}")

    summary = (
        df.groupby("rnaseq_cluster")["mutation_burden"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    print("\nMutation burden summary by cluster:")
    print(summary.to_string(index=False))

    print(f"\nKruskal-Wallis statistic: {statistic:.4f}")
    print(f"Kruskal-Wallis p-value: {p_value:.4e}")


def save_summary(df: pd.DataFrame, statistic: float, p_value: float) -> None:
    out_dir = get_path_from_root("reports/tables/mutations")
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = (
        df.groupby("rnaseq_cluster")["mutation_burden"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )

    summary["kruskal_statistic"] = statistic
    summary["kruskal_p_value"] = p_value

    out_file = out_dir / "mutation_burden_by_rnaseq_cluster.csv"
    summary.to_csv(out_file, index=False)

    print(f"\n[OK] Saved burden summary: {out_file}")


def plot_burden(df: pd.DataFrame, p_value: float) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    clusters = sorted(df["rnaseq_cluster"].unique())
    data = [
        df.loc[df["rnaseq_cluster"] == cluster, "mutation_burden"]
        for cluster in clusters
    ]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.boxplot(data, labels=[f"Cluster {c}" for c in clusters], showfliers=False)

    ax.set_title("Mutation burden by RNA-seq cluster")
    ax.set_xlabel("RNA-seq cluster")
    ax.set_ylabel("Non-silent mutation burden")
    ax.text(
        0.03,
        0.95,
        f"Kruskal-Wallis p = {p_value:.2e}",
        transform=ax.transAxes,
        va="top",
    )

    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_mutation_burden_by_rnaseq_cluster.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved mutation burden figure: {out_file}")


def main() -> int:
    try:
        clusters = load_clusters()
        burden = load_burden()

        df = merge_data(clusters, burden)
        statistic, p_value = run_test(df)

        summarize(df, statistic, p_value)
        save_summary(df, statistic, p_value)
        plot_burden(df, p_value)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())