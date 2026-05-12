from __future__ import annotations

import sys
import pandas as pd

from lifelines import CoxPHFitter

from multiomics_uc.paths import get_path_from_root


def load_data() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    cols = [
        "patient_id",
        "survival_time",
        "survival_event",
        "age_at_index",
        "gender",
        "rnaseq_cluster",
    ]

    df = df[cols].copy()

    df["survival_time"] = pd.to_numeric(df["survival_time"], errors="coerce")
    df["survival_event"] = pd.to_numeric(df["survival_event"], errors="coerce")
    df["age_at_index"] = pd.to_numeric(df["age_at_index"], errors="coerce")

    df = df.dropna()
    df = df[df["survival_time"] > 0].copy()

    df["gender_male"] = (df["gender"].str.lower() == "male").astype(int)

    cluster_dummies = pd.get_dummies(
        df["rnaseq_cluster"].astype(str),
        prefix="cluster",
        drop_first=True,
    )

    model_df = pd.concat(
        [
            df[
                [
                    "patient_id",
                    "survival_time",
                    "survival_event",
                    "age_at_index",
                    "gender_male",
                ]
            ],
            cluster_dummies,
        ],
        axis=1,
    )

    return model_df


def fit_cox_model(df: pd.DataFrame, covariates: list[str]) -> CoxPHFitter:
    model_df = df[["survival_time", "survival_event"] + covariates].copy()

    cph = CoxPHFitter()
    cph.fit(
        model_df,
        duration_col="survival_time",
        event_col="survival_event",
    )

    return cph


def model_summary(cph: CoxPHFitter, model_name: str) -> pd.DataFrame:
    summary = cph.summary.reset_index()
    summary = summary.rename(columns={"covariate": "variable"})
    summary["model"] = model_name
    summary["concordance_index"] = cph.concordance_index_
    summary["partial_AIC"] = cph.AIC_partial_

    return summary


def summarize_results(
    clinical_model: CoxPHFitter,
    cluster_model: CoxPHFitter,
) -> None:
    print("\n Cox Survival Baseline Summary")

    print("\nClinical model:")
    print(f"C-index: {clinical_model.concordance_index_:.3f}")
    print(f"Partial AIC: {clinical_model.AIC_partial_:.2f}")

    print("\nClinical + RNA-seq cluster model:")
    print(f"C-index: {cluster_model.concordance_index_:.3f}")
    print(f"Partial AIC: {cluster_model.AIC_partial_:.2f}")

    delta_c = cluster_model.concordance_index_ - clinical_model.concordance_index_
    delta_aic = cluster_model.AIC_partial_ - clinical_model.AIC_partial_

    print(f"\nDelta C-index: {delta_c:.3f}")
    print(f"Delta partial AIC: {delta_aic:.2f}")


def save_results(results: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/survival_models")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "cox_cluster_baseline_results.csv"
    results.to_csv(out_file, index=False)

    print(f"\n[OK] Saved Cox model results: {out_file}")


def main() -> int:
    try:
        raw = load_data()
        df = prepare_data(raw)

        clinical_covariates = [
            "age_at_index",
            "gender_male",
        ]

        cluster_covariates = [
            col for col in df.columns if col.startswith("cluster_")
        ]

        clinical_plus_cluster_covariates = clinical_covariates + cluster_covariates

        clinical_model = fit_cox_model(df, clinical_covariates)
        cluster_model = fit_cox_model(df, clinical_plus_cluster_covariates)

        summarize_results(clinical_model, cluster_model)

        results = pd.concat(
            [
                model_summary(clinical_model, "clinical_only"),
                model_summary(cluster_model, "clinical_plus_rnaseq_cluster"),
            ],
            axis=0,
            ignore_index=True,
        )

        save_results(results)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())