from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import chi2_contingency, kruskal

from multiomics_uc.paths import get_path_from_root


CATEGORICAL_VARS = [
    "gender",
    "vital_status",
    "ajcc_pathologic_stage",
    "ajcc_pathologic_t",
    "ajcc_pathologic_n",
    "ajcc_pathologic_m",
]

CONTINUOUS_VARS = [
    "age_at_index",
    "survival_time",
]


def load_data() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)


def clean_stage(value):
    if pd.isna(value):
        return pd.NA

    value = str(value).strip()

    stage_map = {
        "Stage 0a": "Stage 0",
        "Stage 0is": "Stage 0",
        "Stage I": "Stage I",
        "Stage II": "Stage II",
        "Stage IIB": "Stage II",
        "Stage IIC": "Stage II",
        "Stage III": "Stage III",
        "Stage IIIA": "Stage III",
        "Stage IV": "Stage IV",
    }

    return stage_map.get(value, value)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["rnaseq_cluster"] = df["rnaseq_cluster"].astype(str)

    if "ajcc_pathologic_stage" in df.columns:
        df["ajcc_pathologic_stage_clean"] = df["ajcc_pathologic_stage"].apply(clean_stage)

    df["age_at_index"] = pd.to_numeric(df["age_at_index"], errors="coerce")
    df["survival_time"] = pd.to_numeric(df["survival_time"], errors="coerce")
    df["survival_event"] = pd.to_numeric(df["survival_event"], errors="coerce")

    return df


def analyze_categorical(df: pd.DataFrame, variable: str) -> dict:
    table = pd.crosstab(df["rnaseq_cluster"], df[variable], dropna=False)

    chi2, p_value, dof, expected = chi2_contingency(table)

    return {
        "variable": variable,
        "test": "chi_square",
        "statistic": chi2,
        "p_value": p_value,
        "degrees_of_freedom": dof,
        "n": int(table.values.sum()),
    }


def analyze_continuous(df: pd.DataFrame, variable: str) -> dict:
    groups = []

    for _, subset in df.groupby("rnaseq_cluster"):
        values = pd.to_numeric(subset[variable], errors="coerce").dropna()
        if len(values) > 0:
            groups.append(values)

    statistic, p_value = kruskal(*groups)

    return {
        "variable": variable,
        "test": "kruskal_wallis",
        "statistic": statistic,
        "p_value": p_value,
        "degrees_of_freedom": len(groups) - 1,
        "n": int(df[variable].notna().sum()),
    }


def save_association_results(results: list[dict]) -> None:
    out_dir = get_path_from_root("reports/tables/clinical_associations")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "rnaseq_cluster_clinical_association_tests.csv"
    pd.DataFrame(results).to_csv(out_file, index=False)

    print(f"[OK] Saved association tests: {out_file}")


def save_contingency_tables(df: pd.DataFrame, categorical_vars: list[str]) -> None:
    out_dir = get_path_from_root("reports/tables/clinical_associations")
    out_dir.mkdir(parents=True, exist_ok=True)

    for var in categorical_vars:
        if var not in df.columns:
            continue

        table = pd.crosstab(df["rnaseq_cluster"], df[var], dropna=False)
        out_file = out_dir / f"contingency_{var}.csv"
        table.to_csv(out_file)

        print(f"[OK] Saved contingency table for {var}: {out_file}")


def plot_stage_distribution(df: pd.DataFrame) -> None:
    if "ajcc_pathologic_stage_clean" not in df.columns:
        return

    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    table = pd.crosstab(
        df["rnaseq_cluster"],
        df["ajcc_pathologic_stage_clean"],
        normalize="index",
        dropna=False,
    )

    ax = table.plot(
        kind="bar",
        stacked=True,
        figsize=(9, 6),
    )

    ax.set_title("Pathologic stage distribution by RNA-seq cluster")
    ax.set_xlabel("RNA-seq cluster")
    ax.set_ylabel("Proportion of patients")
    ax.legend(title="Stage", bbox_to_anchor=(1.05, 1), loc="upper left")

    plt.tight_layout()

    out_file = figures_dir / "tcga_blca_stage_distribution_by_rnaseq_cluster.png"
    plt.savefig(out_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"[OK] Saved stage distribution figure: {out_file}")


def plot_age_distribution(df: pd.DataFrame) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    clusters = sorted(df["rnaseq_cluster"].dropna().unique())
    data = [
        df.loc[df["rnaseq_cluster"] == cluster, "age_at_index"].dropna()
        for cluster in clusters
    ]

    ax.boxplot(data, labels=[f"Cluster {c}" for c in clusters])
    ax.set_title("Age distribution by RNA-seq cluster")
    ax.set_xlabel("RNA-seq cluster")
    ax.set_ylabel("Age at index")

    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_age_distribution_by_rnaseq_cluster.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved age distribution figure: {out_file}")


def summarize(df: pd.DataFrame, results: list[dict]) -> None:
    print("\n Clinical Association Summary")
    print(f"Patients: {df['patient_id'].nunique()}")

    print("\nCluster sizes:")
    print(df["rnaseq_cluster"].value_counts().sort_index())

    print("\nAssociation tests:")
    res_df = pd.DataFrame(results).sort_values("p_value")
    print(res_df.to_string(index=False))


def main() -> int:
    try:
        df = load_data()
        df = prepare_data(df)

        categorical_vars = [
            "gender",
            "vital_status",
            "ajcc_pathologic_stage_clean",
            "ajcc_pathologic_t",
            "ajcc_pathologic_n",
            "ajcc_pathologic_m",
        ]

        results = []

        for var in categorical_vars:
            if var in df.columns:
                test_df = df.dropna(subset=[var])
                if test_df[var].nunique() > 1:
                    results.append(analyze_categorical(test_df, var))

        for var in CONTINUOUS_VARS:
            if var in df.columns:
                results.append(analyze_continuous(df, var))

        summarize(df, results)
        save_association_results(results)
        save_contingency_tables(df, categorical_vars)
        plot_stage_distribution(df)
        plot_age_distribution(df)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())