from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from lifelines import KaplanMeierFitter
from lifelines.statistics import multivariate_logrank_test

from multiomics_uc.paths import get_path_from_root


def load_clustered_data() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)


def prepare_survival_data(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = ["patient_id", "survival_time", "survival_event", "rnaseq_cluster"]

    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")

    survival_df = df[required_cols].copy()

    survival_df["survival_time"] = pd.to_numeric(
        survival_df["survival_time"],
        errors="coerce",
    )

    survival_df["survival_event"] = pd.to_numeric(
        survival_df["survival_event"],
        errors="coerce",
    )

    survival_df = survival_df.dropna(
        subset=["survival_time", "survival_event", "rnaseq_cluster"]
    )

    survival_df = survival_df[survival_df["survival_time"] > 0]

    return survival_df


def run_logrank_test(survival_df: pd.DataFrame):
    result = multivariate_logrank_test(
        event_durations=survival_df["survival_time"],
        groups=survival_df["rnaseq_cluster"],
        event_observed=survival_df["survival_event"],
    )

    return result


def plot_kaplan_meier(survival_df: pd.DataFrame, p_value: float) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))

    kmf = KaplanMeierFitter()

    for cluster, subset in survival_df.groupby("rnaseq_cluster"):
        kmf.fit(
            durations=subset["survival_time"],
            event_observed=subset["survival_event"],
            label=f"Cluster {cluster} (n={len(subset)})",
        )
        kmf.plot_survival_function(ax=ax, ci_show=True)

    ax.set_title("Overall survival by RNA-seq transcriptomic cluster")
    ax.set_xlabel("Time, days")
    ax.set_ylabel("Survival probability")
    ax.text(
        0.02,
        0.05,
        f"Log-rank p = {p_value:.3e}",
        transform=ax.transAxes,
        fontsize=10,
    )
    ax.legend(frameon=False)
    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_kaplan_meier_by_rnaseq_cluster.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved Kaplan-Meier figure: {out_file}")


def save_summary(survival_df: pd.DataFrame, p_value: float, test_statistic: float) -> None:
    reports_dir = get_path_from_root("reports/tables")
    reports_dir.mkdir(parents=True, exist_ok=True)

    cluster_summary = (
        survival_df.groupby("rnaseq_cluster")
        .agg(
            n_patients=("patient_id", "count"),
            n_events=("survival_event", "sum"),
            median_survival_time=("survival_time", "median"),
            mean_survival_time=("survival_time", "mean"),
        )
        .reset_index()
    )

    cluster_summary["global_logrank_p_value"] = p_value
    cluster_summary["global_logrank_test_statistic"] = test_statistic

    out_file = reports_dir / "tcga_blca_survival_by_rnaseq_cluster_summary.csv"
    cluster_summary.to_csv(out_file, index=False)

    print(f"[OK] Saved survival summary: {out_file}")


def summarize(survival_df: pd.DataFrame, p_value: float, test_statistic: float) -> None:
    print("\n Survival by RNA-seq Cluster Summary")
    print(f"Patients included: {len(survival_df)}")
    print("\nCluster sizes:")
    print(survival_df["rnaseq_cluster"].value_counts().sort_index())

    print("\nEvents by cluster:")
    print(survival_df.groupby("rnaseq_cluster")["survival_event"].sum())

    print(f"\nGlobal log-rank test statistic: {test_statistic:.4f}")
    print(f"Global log-rank p-value: {p_value:.4e}")


def main() -> int:
    try:
        clustered = load_clustered_data()
        survival_df = prepare_survival_data(clustered)

        logrank_result = run_logrank_test(survival_df)

        p_value = float(logrank_result.p_value)
        test_statistic = float(logrank_result.test_statistic)

        summarize(survival_df, p_value, test_statistic)
        save_summary(survival_df, p_value, test_statistic)
        plot_kaplan_meier(survival_df, p_value)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())