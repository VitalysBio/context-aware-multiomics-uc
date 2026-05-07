from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from multiomics_uc.paths import get_path_from_root


def find_first_star_counts_file() -> Path:
    raw_dir = get_path_from_root("data/raw/tcga_blca/rna_seq")

    files = list(raw_dir.glob("*/*.tsv"))

    if not files:
        raise FileNotFoundError(f"No STAR count files found under: {raw_dir}")

    return files[0]


def load_gene_annotation(file_path: Path) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t", skiprows=1)

    required_cols = {"gene_id", "gene_name", "gene_type"}
    missing = required_cols - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns {missing} in file: {file_path}")

    annotation = df[["gene_id", "gene_name", "gene_type"]].copy()

    annotation = annotation[
        ~annotation["gene_id"].astype(str).str.startswith("N_")
    ].copy()

    annotation["gene_id_clean"] = (
        annotation["gene_id"]
        .astype(str)
        .str.replace(r"\.\d+$", "", regex=True)
    )

    annotation = annotation.drop_duplicates(subset=["gene_id"])

    return annotation


def summarize(annotation: pd.DataFrame, source_file: Path) -> None:
    print("\n Gene Annotation Extraction Summary")
    print(f"Source file: {source_file}")
    print(f"Genes annotated: {len(annotation)}")
    print(f"Unique gene IDs: {annotation['gene_id'].nunique()}")
    print(f"Unique gene names: {annotation['gene_name'].nunique()}")

    print("\nTop gene types:")
    print(annotation["gene_type"].value_counts().head(15))


def save(annotation: pd.DataFrame) -> None:
    out_dir = get_path_from_root("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "tcga_blca_gene_annotation_table.csv"
    annotation.to_csv(out_file, index=False)

    print(f"\n[OK] Saved gene annotation table: {out_file}")


def main() -> int:
    try:
        source_file = find_first_star_counts_file()
        annotation = load_gene_annotation(source_file)

        summarize(annotation, source_file)
        save(annotation)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())