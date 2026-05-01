from __future__ import annotations

import sys
import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_master_table() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_master_clinical_table.csv")
    return pd.read_csv(path)


def summarize_missingness(df: pd.DataFrame) -> pd.DataFrame:
    rows = []

    for col in df.columns:
        missing = df[col].isna().sum()
        rows.append(
            {
                "column": col,
                "missing_count": missing,
                "missing_percent": round(100 * missing / len(df), 2),
                "non_missing_count": len(df) - missing,
            }
        )

    return pd.DataFrame(rows).sort_values("missing_percent", ascending=False)


def summarize_key_variables(df: pd.DataFrame) -> None:
    print("\n Clinical Audit Summary ")
    print(f"Rows: {len(df)}")
    print(f"Unique patients: {df['patient_id'].nunique()}")

    key_cols = [
        "gender",
        "vital_status",
        "survival_time",
        "survival_event",
        "age_at_index",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
        "primary_diagnosis",
    ]

    print("\nMissingness in key variables:")
    for col in key_cols:
        if col in df.columns:
            missing = df[col].isna().sum()
            print(f"{col}: {missing} missing ({100 * missing / len(df):.1f}%)")

    print("\nVital status distribution:")
    print(df["vital_status"].value_counts(dropna=False))

    print("\nPathologic stage distribution:")
    if "ajcc_pathologic_stage" in df.columns:
        print(df["ajcc_pathologic_stage"].value_counts(dropna=False))

    usable_survival = df["survival_time"].notna() & df["survival_event"].notna()
    print("\nPatients usable for survival analysis:")
    print(f"{usable_survival.sum()} / {len(df)}")


def save_missingness_report(report: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "tcga_blca_clinical_missingness_report.csv"
    report.to_csv(out_file, index=False)

    print(f"\n[OK] Saved: {out_file}")


def main() -> int:
    try:
        df = load_master_table()
        summarize_key_variables(df)

        report = summarize_missingness(df)
        save_missingness_report(report)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())