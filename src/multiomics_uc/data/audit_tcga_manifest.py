from __future__ import annotations

import pandas as pd
from pathlib import Path
import sys

from multiomics_uc.paths import get_path_from_root


def load_manifest(manifest_path: Path) -> pd.DataFrame:
    if not manifest_path.exists():
        raise FileNotFoundError(f"Manifest not found: {manifest_path}")

    df = pd.read_csv(manifest_path, sep="\t")
    return df


def validate_manifest(df: pd.DataFrame) -> None:
    required_columns = {"id", "filename", "md5", "size", "state"}

    missing = required_columns - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def summarize_manifest(df: pd.DataFrame) -> None:
    print("\n=== Manifest Summary ===")
    print(f"Total files: {len(df)}")

    total_size_gb = df["size"].sum() / (1024**3)
    print(f"Total size (GB): {total_size_gb:.2f}")

    print("\nState distribution:")
    print(df["state"].value_counts())


def extract_sample_ids(df: pd.DataFrame) -> pd.Series:
    # TCGA filenames contain sample barcodes
    # example: XXXXXXX.rna_seq.augmented_star_gene_counts.tsv

    sample_ids = df["filename"].str.split(".").str[0]

    print("\nUnique sample IDs:", sample_ids.nunique())

    return sample_ids


def save_report(df: pd.DataFrame, sample_ids: pd.Series) -> None:
    output_dir = get_path_from_root("reports", "tables")
    output_dir.mkdir(parents=True, exist_ok=True)

    df_out = df.copy()
    df_out["sample_id"] = sample_ids

    output_file = output_dir / "tcga_blca_rna_seq_manifest_audit.csv"
    df_out.to_csv(output_file, index=False)

    print(f"\n[OK] Report saved to: {output_file}")


def main() -> int:
    try:
        manifest_path = get_path_from_root(
            "data",
            "external",
            "tcga_blca",
            "manifests",
            "rna_seq",
            "gdc_manifest_tcga_blca_rna_seq_star_counts.txt",
        )

        df = load_manifest(manifest_path)

        validate_manifest(df)
        summarize_manifest(df)

        sample_ids = extract_sample_ids(df)

        save_report(df, sample_ids)

        print("\n[OK] Manifest audit completed successfully")
        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())