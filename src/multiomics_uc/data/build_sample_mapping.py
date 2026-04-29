from __future__ import annotations

import pandas as pd
from pathlib import Path
import json
import sys

from multiomics_uc.paths import get_path_from_root


def load_manifest():
    path = get_path_from_root(
        "data/external/tcga_blca/manifests/rna_seq",
        "gdc_manifest_tcga_blca_rna_seq_star_counts.txt"
    )
    return pd.read_csv(path, sep="\t")


def load_sample_sheet():
    path = get_path_from_root(
        "data/external/tcga_blca/metadata",
        "sample_sheet.tsv"
    )
    return pd.read_csv(path, sep="\t")


def extract_patient_id(sample_barcode: str) -> str:
    # TCGA barcode structure:
    # TCGA-XX-XXXX-...
    return "-".join(sample_barcode.split("-")[0:3])


def build_mapping():
    manifest = load_manifest()
    sample_sheet = load_sample_sheet()

    # rename for clarity
    sample_sheet = sample_sheet.rename(columns={
        "File Name": "filename",
        "Sample ID": "sample_barcode"
    })

    df = manifest.merge(sample_sheet, on="filename", how="left")

    df["patient_id"] = df["sample_barcode"].apply(
        lambda x: extract_patient_id(x) if pd.notnull(x) else None
    )

    return df


def summarize(df: pd.DataFrame):
    print("\n Mapping Summary ")
    print(f"Files: {len(df)}")
    print(f"Unique samples: {df['sample_barcode'].nunique()}")
    print(f"Unique patients: {df['patient_id'].nunique()}")

    print("\nSamples per patient (top 10):")
    print(df.groupby("patient_id")["sample_barcode"].nunique().sort_values(ascending=False).head(10))


def save(df: pd.DataFrame):
    out = get_path_from_root("reports/tables")
    out.mkdir(parents=True, exist_ok=True)

    file = out / "tcga_blca_sample_mapping.csv"
    df.to_csv(file, index=False)

    print(f"\n[OK] Saved: {file}")


def main():
    try:
        df = build_mapping()
        summarize(df)
        save(df)
        return 0
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())