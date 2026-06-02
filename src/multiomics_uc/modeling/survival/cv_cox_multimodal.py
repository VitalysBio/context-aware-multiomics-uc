from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from lifelines import CoxPHFitter
from lifelines.utils import concordance_index
from sklearn.model_selection import KFold

from multiomics_uc.paths import get_path_from_root


FEATURES = [
    "age_at_index",
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PC6", "PC7", "PC8", "PC9", "PC10",
    "mutation_burden",
    "TP53_mut", "FGFR3_mut", "RB1_mut",
    "KDM6A_mut", "ARID1A_mut", "PIK3CA_mut",
    "EGFR_amp", "ERBB2_amp", "MYC_amp",
    "CCND1_amp", "CDKN2A_del", "RB1_del",
]


def load_data() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/integrated_patient_feature_table_full_multiomics.csv"
    )
    return pd.read_csv(path)


def prepare_data(df: pd.DataFrame) -> pd.DataFrame:
    cols = ["patient_id", "survival_time", "survival_event"] + FEATURES
    df = df[cols].copy()

    for col in ["survival_time", "survival_event"] + FEATURES:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna()
    df = df[df["survival_time"] > 0].copy()

    return df


def fit_cox(train_df: pd.DataFrame) -> CoxPHFitter:
    model_df = train_df[["survival_time", "survival_event"] + FEATURES].copy()

    cph = CoxPHFitter(penalizer=0.1)
    cph.fit(
        model_df,
        duration_col="survival_time",
        event_col="survival_event",
    )

    return cph


def run_cross_validation(df: pd.DataFrame, n_splits: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    kf = KFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=42,
    )

    fold_rows = []
    prediction_rows = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(df), start=1):
        train_df = df.iloc[train_idx].copy()
        test_df = df.iloc[test_idx].copy()

        cph = fit_cox(train_df)

        partial_hazard = cph.predict_partial_hazard(
            test_df[FEATURES]
        ).values.ravel()

        c_index = concordance_index(
            event_times=test_df["survival_time"],
            predicted_scores=-partial_hazard,
            event_observed=test_df["survival_event"],
        )

        fold_rows.append(
            {
                "fold": fold,
                "n_train": len(train_df),
                "n_test": len(test_df),
                "c_index": float(c_index),
            }
        )

        for _, row in test_df.iterrows():
            prediction_rows.append(
                {
                    "patient_id": row["patient_id"],
                    "fold": fold,
                    "risk_score": float(
                        partial_hazard[test_df.index.get_loc(row.name)]
                    ),
                    "survival_time": row["survival_time"],
                    "survival_event": row["survival_event"],
                }
            )

        print(f"Fold {fold}: C-index = {c_index:.3f}")

    return pd.DataFrame(fold_rows), pd.DataFrame(prediction_rows)


def summarize(fold_results: pd.DataFrame) -> None:
    print("\nCross-validated Penalized Cox Full Multiomics")
    print(fold_results.to_string(index=False))

    print("\nMean C-index:")
    print(fold_results["c_index"].mean().round(3))

    print("\nStd C-index:")
    print(fold_results["c_index"].std().round(3))


def save_outputs(fold_results: pd.DataFrame, predictions: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/survival_models")
    out_dir.mkdir(parents=True, exist_ok=True)

    fold_file = out_dir / "cox_full_multiomics_cv_results.csv"
    pred_file = out_dir / "cox_full_multiomics_cv_predictions.csv"

    fold_results.to_csv(fold_file, index=False)
    predictions.to_csv(pred_file, index=False)

    print(f"\n[OK] Saved Cox CV results: {fold_file}")
    print(f"[OK] Saved Cox CV predictions: {pred_file}")


def main() -> int:
    try:
        raw = load_data()
        df = prepare_data(raw)

        print("\nCross-validated Cox Full Multiomics")
        print(f"Patients: {df.shape[0]}")
        print(f"Features: {len(FEATURES)}")

        fold_results, predictions = run_cross_validation(df, n_splits=5)

        summarize(fold_results)
        save_outputs(fold_results, predictions)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())