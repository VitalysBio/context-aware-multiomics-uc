from __future__ import annotations

import sys
import pandas as pd
from lifelines import CoxPHFitter

from multiomics_uc.paths import get_path_from_root


CLINICAL = ["age_at_index"]
RNA = ["rnaseq_cluster"] + [f"PC{i}" for i in range(1, 11)]
MUTATION = [
    "mutation_burden",
    "TP53_mut",
    "FGFR3_mut",
    "RB1_mut",
    "KDM6A_mut",
    "ARID1A_mut",
    "PIK3CA_mut",
]
CNV = [
    "EGFR_amp",
    "ERBB2_amp",
    "MYC_amp",
    "CCND1_amp",
    "CDKN2A_del",
    "RB1_del",
]


MODEL_SPECS = {
    "clinical_only": CLINICAL,
    "rna_only": RNA,
    "mutation_only": MUTATION,
    "clinical_plus_rna": CLINICAL + RNA,
    "clinical_plus_mutation": CLINICAL + MUTATION,
    "clinical_plus_rna_mutation": CLINICAL + RNA + MUTATION,
    "full_multiomics": CLINICAL + RNA + MUTATION + CNV,
}


def load_clean_table() -> pd.DataFrame:
    path = get_path_from_root("data/processed/integrated_patient_feature_table_clean.csv")
    return pd.read_csv(path)


def prepare_model_data(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    required = ["survival_time", "survival_event"] + features
    model_df = df[[c for c in required if c in df.columns]].copy()

    if "rnaseq_cluster" in model_df.columns:
        model_df["rnaseq_cluster"] = model_df["rnaseq_cluster"].astype(str)
        model_df = pd.get_dummies(
            model_df,
            columns=["rnaseq_cluster"],
            prefix="cluster",
            drop_first=True,
        )

    for col in model_df.columns:
        if col not in ["survival_time", "survival_event"]:
            model_df[col] = pd.to_numeric(model_df[col], errors="coerce")

    model_df = model_df.dropna()
    model_df = model_df[model_df["survival_time"] > 0]

    return model_df


def fit_cox(model_df: pd.DataFrame) -> CoxPHFitter:
    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(
        model_df,
        duration_col="survival_time",
        event_col="survival_event",
    )
    return cph


def run_benchmark(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_rows = []
    coefficient_tables = []

    for model_name, features in MODEL_SPECS.items():
        model_df = prepare_model_data(df, features)

        if model_df.shape[0] < 50:
            print(f"[SKIP] {model_name}: too few patients after dropna")
            continue

        cph = fit_cox(model_df)

        summary_rows.append(
            {
                "model": model_name,
                "n_patients": model_df.shape[0],
                "n_features": model_df.shape[1] - 2,
                "c_index": cph.concordance_index_,
                "partial_AIC": cph.AIC_partial_,
            }
        )

        coef = cph.summary.reset_index().rename(columns={"covariate": "feature"})
        coef["model"] = model_name
        coefficient_tables.append(coef)

    benchmark = pd.DataFrame(summary_rows).sort_values("c_index", ascending=False)
    coefficients = pd.concat(coefficient_tables, axis=0, ignore_index=True)

    return benchmark, coefficients


def summarize(benchmark: pd.DataFrame) -> None:
    print("\nCox Multimodal Survival Benchmark")
    print(benchmark.to_string(index=False))


def save_outputs(benchmark: pd.DataFrame, coefficients: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/survival_models")
    out_dir.mkdir(parents=True, exist_ok=True)

    benchmark_file = out_dir / "cox_multimodal_benchmark_summary.csv"
    coefficients_file = out_dir / "cox_multimodal_benchmark_coefficients.csv"

    benchmark.to_csv(benchmark_file, index=False)
    coefficients.to_csv(coefficients_file, index=False)

    print(f"\n[OK] Saved benchmark summary: {benchmark_file}")
    print(f"[OK] Saved coefficients: {coefficients_file}")


def main() -> int:
    try:
        df = load_clean_table()

        benchmark, coefficients = run_benchmark(df)

        summarize(benchmark)
        save_outputs(benchmark, coefficients)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())