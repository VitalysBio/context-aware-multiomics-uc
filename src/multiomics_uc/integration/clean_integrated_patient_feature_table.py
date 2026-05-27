from __future__ import annotations

import sys
import pandas as pd

from multiomics_uc.paths import get_path_from_root


CNV_COLUMNS = [
    "EGFR_amp",
    "ERBB2_amp",
    "MYC_amp",
    "CCND1_amp",
    "CDKN2A_del",
    "RB1_del",
]


def load_integrated_table() -> pd.DataFrame:
    path = get_path_from_root("data/processed/integrated_patient_feature_table.csv")
    return pd.read_csv(path)


def clean_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}

    preferred_columns = {
        "age_at_index_x": "age_at_index",
        "gender_x": "gender",
        "ajcc_pathologic_stage_x": "ajcc_pathologic_stage",
        "survival_time_x": "survival_time",
        "survival_event_x": "survival_event",
    }

    for old, new in preferred_columns.items():
        if old in df.columns:
            rename_map[old] = new

    df = df.rename(columns=rename_map)

    drop_cols = [
        "age_at_index_y",
        "gender_y",
        "ajcc_pathologic_stage_y",
        "survival_time_y",
        "survival_event_y",
    ]

    drop_cols = [col for col in drop_cols if col in df.columns]
    df = df.drop(columns=drop_cols)

    return df


def add_availability_flags(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    mutation_cols = [col for col in df.columns if col.endswith("_mut")]

    df["has_mutation_data"] = df[mutation_cols].notna().any(axis=1)

    existing_cnv_cols = [col for col in CNV_COLUMNS if col in df.columns]
    df["has_cnv_data"] = df[existing_cnv_cols].notna().any(axis=1)

    df["has_full_multiomics"] = (
        df["has_mutation_data"] & df["has_cnv_data"]
    )

    return df


def fill_mutation_missing(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    mutation_cols = [col for col in df.columns if col.endswith("_mut")]

    for col in mutation_cols:
        df[col] = df[col].fillna(0).astype(int)
    
    if "mutation_burden" in df.columns:
        df["has_mutation_burden"] = (df["mutation_burden"].notna())

    return df


def build_full_multiomics_subset(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["has_full_multiomics"]].copy()


def summarize(clean_df: pd.DataFrame, full_df: pd.DataFrame) -> None:
    print("\n Clean Integrated Feature Table Summary")
    print(f"Clean table shape: {clean_df.shape}")
    print(f"Patients: {clean_df['patient_id'].nunique()}")

    print("\nData availability:")
    print(clean_df[["has_mutation_data", "has_cnv_data", "has_full_multiomics"]].sum())

    print("\nFull multi-omics table:")
    print(f"Shape: {full_df.shape}")
    print(f"Patients: {full_df['patient_id'].nunique()}")

    print("\nTop missingness after cleaning:")
    print(clean_df.isna().mean().sort_values(ascending=False).head(15))


def save_outputs(clean_df: pd.DataFrame, full_df: pd.DataFrame) -> None:
    processed_dir = get_path_from_root("data/processed")

    clean_file = processed_dir / "integrated_patient_feature_table_clean.csv"
    full_file = processed_dir / "integrated_patient_feature_table_full_multiomics.csv"

    clean_df.to_csv(clean_file, index=False)
    full_df.to_csv(full_file, index=False)

    print(f"\n[OK] Saved clean table: {clean_file}")
    print(f"[OK] Saved full multi-omics table: {full_file}")


def main() -> int:
    try:
        df = load_integrated_table()

        clean_df = clean_duplicate_columns(df)
        clean_df = add_availability_flags(clean_df)
        clean_df = fill_mutation_missing(clean_df)

        full_df = build_full_multiomics_subset(clean_df)

        summarize(clean_df, full_df)
        save_outputs(clean_df, full_df)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())