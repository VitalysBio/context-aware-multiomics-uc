from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sklearn.inspection import permutation_importance

from multiomics_uc.paths import get_path_from_root


FEATURES = [
    "age_at_index",
    "PC1",
    "PC2",
    "PC3",
    "PC4",
    "PC5",
    "PC6",
    "PC7",
    "PC8",
    "PC9",
    "PC10",
    "mutation_burden",
    "TP53_mut",
    "FGFR3_mut",
    "RB1_mut",
    "KDM6A_mut",
    "ARID1A_mut",
    "PIK3CA_mut",
    "EGFR_amp",
    "ERBB2_amp",
    "MYC_amp",
    "CCND1_amp",
    "CDKN2A_del",
    "RB1_del",
]


def load_data() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/integrated_patient_feature_table_full_multiomics.csv"
    )
    return pd.read_csv(path)


def prepare_data(df: pd.DataFrame) -> tuple[pd.DataFrame, object]:
    required = ["survival_time", "survival_event"] + FEATURES

    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    model_df = df[required].copy()

    for col in FEATURES:
        model_df[col] = pd.to_numeric(model_df[col], errors="coerce")

    model_df["survival_time"] = pd.to_numeric(
        model_df["survival_time"],
        errors="coerce",
    )

    model_df["survival_event"] = pd.to_numeric(
        model_df["survival_event"],
        errors="coerce",
    )

    model_df = model_df.dropna()
    model_df = model_df[model_df["survival_time"] > 0].copy()

    X = model_df[FEATURES]

    y = model_df[["survival_event", "survival_time"]].copy()
    y["survival_event"] = y["survival_event"].astype(bool)

    y_structured = y.to_records(index=False)

    return X, y_structured


def fit_model(X: pd.DataFrame, y) -> RandomSurvivalForest:
    model = RandomSurvivalForest(
        n_estimators=1000,
        min_samples_split=10,
        min_samples_leaf=5,
        max_features="sqrt",
        random_state=42,
        n_jobs=-1,
    )

    model.fit(X, y)

    return model


def compute_c_index(model: RandomSurvivalForest, X: pd.DataFrame, y) -> float:
    risk_scores = model.predict(X)

    c_index = concordance_index_censored(
        y["survival_event"],
        y["survival_time"],
        risk_scores,
    )[0]

    return float(c_index)


def compute_permutation_importance(
    model: RandomSurvivalForest,
    X: pd.DataFrame,
    y,
) -> pd.DataFrame:
    result = permutation_importance(
        model,
        X,
        y,
        n_repeats=20,
        random_state=42,
        n_jobs=-1,
    )

    importance = pd.DataFrame(
        {
            "feature": X.columns,
            "importance_mean": result.importances_mean,
            "importance_std": result.importances_std,
        }
    ).sort_values("importance_mean", ascending=False)

    return importance


def save_outputs(c_index: float, importance: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/survival_models")
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame(
        [
            {
                "model": "random_survival_forest_full_multiomics",
                "n_features": len(FEATURES),
                "c_index_apparent": c_index,
            }
        ]
    )

    summary_file = out_dir / "rsf_full_multiomics_summary.csv"
    importance_file = out_dir / "rsf_full_multiomics_permutation_importance.csv"

    summary.to_csv(summary_file, index=False)
    importance.to_csv(importance_file, index=False)

    print(f"\n[OK] Saved RSF summary: {summary_file}")
    print(f"[OK] Saved RSF feature importance: {importance_file}")


def plot_importance(importance: pd.DataFrame) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    top = importance.head(20).sort_values("importance_mean", ascending=True)

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.barh(
        top["feature"],
        top["importance_mean"],
        xerr=top["importance_std"],
    )

    ax.set_title("Random Survival Forest feature importance")
    ax.set_xlabel("Permutation importance")
    ax.set_ylabel("Feature")

    fig.tight_layout()

    out_file = figures_dir / "rsf_full_multiomics_feature_importance.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved RSF importance figure: {out_file}")


def summarize(c_index: float, importance: pd.DataFrame, X: pd.DataFrame) -> None:
    print("\nRandom Survival Forest Full Multiomics")
    print(f"Patients used: {X.shape[0]}")
    print(f"Features used: {X.shape[1]}")
    print(f"Apparent C-index: {c_index:.3f}")

    print("\nTop 20 permutation importance features:")
    print(
        importance.head(20).to_string(
            index=False,
        )
    )


def main() -> int:
    try:
        df = load_data()
        X, y = prepare_data(df)

        model = fit_model(X, y)
        c_index = compute_c_index(model, X, y)

        importance = compute_permutation_importance(
            model,
            X,
            y,
        )

        summarize(c_index, importance, X)
        save_outputs(c_index, importance)
        plot_importance(importance)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())