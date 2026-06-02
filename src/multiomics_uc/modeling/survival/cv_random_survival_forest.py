from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from sklearn.model_selection import KFold
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored

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


def prepare_data(df: pd.DataFrame):
    cols = ["patient_id", "survival_time", "survival_event"] + FEATURES
    df = df[cols].copy()

    for col in ["survival_time", "survival_event"] + FEATURES:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df = df.dropna()
    df = df[df["survival_time"] > 0].copy()

    X = df[FEATURES]

    y_df = df[["survival_event", "survival_time"]].copy()
    y_df["survival_event"] = y_df["survival_event"].astype(bool)
    y = y_df.to_records(index=False)

    return df, X, y


def fit_rsf() -> RandomSurvivalForest:
    return RandomSurvivalForest(
        n_estimators=1000,
        min_samples_split=10,
        min_samples_leaf=5,
        max_features="sqrt",
        random_state=42,
        n_jobs=-1,
    )


def run_cross_validation(X: pd.DataFrame, y, n_splits: int = 5) -> tuple[pd.DataFrame, pd.DataFrame]:
    kf = KFold(
        n_splits=n_splits,
        shuffle=True,
        random_state=42,
    )

    fold_rows = []
    prediction_rows = []

    for fold, (train_idx, test_idx) in enumerate(kf.split(X), start=1):
        X_train = X.iloc[train_idx]
        X_test = X.iloc[test_idx]

        y_train = y[train_idx]
        y_test = y[test_idx]

        model = fit_rsf()
        model.fit(X_train, y_train)

        risk_scores = model.predict(X_test)

        c_index = concordance_index_censored(
            y_test["survival_event"],
            y_test["survival_time"],
            risk_scores,
        )[0]

        fold_rows.append(
            {
                "fold": fold,
                "n_train": len(train_idx),
                "n_test": len(test_idx),
                "c_index": float(c_index),
            }
        )

        for idx, risk in zip(test_idx, risk_scores):
            prediction_rows.append(
                {
                    "row_index": int(idx),
                    "fold": fold,
                    "risk_score": float(risk),
                }
            )

        print(f"Fold {fold}: C-index = {c_index:.3f}")

    return pd.DataFrame(fold_rows), pd.DataFrame(prediction_rows)


def summarize(fold_results: pd.DataFrame) -> None:
    print("\n Cross-validated RSF Summary")
    print(fold_results.to_string(index=False))

    print("\nMean C-index:")
    print(fold_results["c_index"].mean().round(3))

    print("\nStd C-index:")
    print(fold_results["c_index"].std().round(3))


def save_outputs(fold_results: pd.DataFrame, predictions: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/survival_models")
    out_dir.mkdir(parents=True, exist_ok=True)

    fold_file = out_dir / "rsf_full_multiomics_cv_results.csv"
    pred_file = out_dir / "rsf_full_multiomics_cv_predictions.csv"

    fold_results.to_csv(fold_file, index=False)
    predictions.to_csv(pred_file, index=False)

    print(f"\n[OK] Saved CV results: {fold_file}")
    print(f"[OK] Saved CV predictions: {pred_file}")


def main() -> int:
    try:
        df = load_data()
        _, X, y = prepare_data(df)

        print("\n Cross-validated Random Survival Forest")
        print(f"Patients: {X.shape[0]}")
        print(f"Features: {X.shape[1]}")

        fold_results, predictions = run_cross_validation(X, y, n_splits=5)

        summarize(fold_results)
        save_outputs(fold_results, predictions)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())