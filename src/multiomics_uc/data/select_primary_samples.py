from __future__ import annotations

import pandas as pd
from pathlib import Path
import sys

from multiomics_uc.paths import get_path_from_root


def load_mapping():
    path = get_path_from_root("reports/tables/tcga_blca_sample_mapping.csv")
    return pd.read_csv(path)


def extract_sample_type(sample_barcode: str) -> str:
    # TCGA barcode format:
    # TCGA-XX-XXXX-YYA
    # YY = sample type
    return sample_barcode.split("-")[3][:2]


def filter_primary(df: pd.DataFrame) -> pd.DataFrame:
    df["sample_type"] = df["sample_barcode"].apply(
        lambda x: extract_sample_type(x) if pd.notnull(x) else None
    )

    primary_df = df[df["sample_type"] == "01"]

    return primary_df


def select_one_per_patient(df: pd.DataFrame) -> pd.DataFrame:
    # in case multiple primary samples exist, keep first deterministically
    df_sorted = df.sort_values(by=["patient_id", "sample_barcode"])

    selected = df_sorted.groupby("patient_id").first().reset_index()

    return selected


def summarize(original: pd.DataFrame, selected: pd.DataFrame):
    print("\n Selection Summary ")
    print(f"Original files: {len(original)}")
    print(f"Original samples: {original['sample_barcode'].nunique()}")
    print(f"Original patients: {original['patient_id'].nunique()}")

    print("\nSample type distribution:")
    print(original["sample_type"].value_counts(dropna=False))

    print(f"\nFinal selected files: {len(selected)}")
    print(f"Final selected samples: {selected['sample_barcode'].nunique()}")
    print(f"Final selected patients: {selected['patient_id'].nunique()}")


def save(df: pd.DataFrame):
    out = get_path_from_root("data/processed")
    out.mkdir(parents=True, exist_ok=True)

    file = out / "tcga_blca_primary_samples.csv"
    df.to_csv(file, index=False)

    print(f"\n[OK] Saved: {file}")


def main():
    try:
        df = load_mapping()

        primary_df = filter_primary(df)
        selected_df = select_one_per_patient(primary_df)

        summarize(df, selected_df)
        save(selected_df)

        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())