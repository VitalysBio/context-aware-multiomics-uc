from __future__ import annotations

import pandas as pd
from pathlib import Path
import sys

from multiomics_uc.paths import get_path_from_root


def load_primary_samples():
    path = get_path_from_root("data/processed/tcga_blca_primary_samples.csv")
    return pd.read_csv(path)


def get_file_path_map():
    raw_dir = get_path_from_root("data/raw/tcga_blca/rna_seq")

    file_map = {}

    for folder in raw_dir.iterdir():
        if folder.is_dir():
            files = list(folder.glob("*.tsv"))
            if len(files) == 1:
                file_map[folder.name] = files[0]

    return file_map


def load_single_file(file_path: Path) -> pd.DataFrame:
    """
    Load one GDC STAR Counts file.

    GDC STAR count files contain a first metadata line:
    gene-model: GENCODE v36

    The actual table starts on the second line.
    """
    df = pd.read_csv(file_path, sep="\t", skiprows=1)

    required_cols = {"gene_id", "unstranded"}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns {missing} in file: {file_path}")

    # Remove technical summary rows such as N_unmapped, N_multimapping, etc.
    df = df[~df["gene_id"].astype(str).str.startswith("N_")]

    # Keep raw unstranded counts for the first expression matrix
    df = df[["gene_id", "unstranded"]].copy()

    # Ensure counts are numeric
    df["unstranded"] = pd.to_numeric(df["unstranded"], errors="coerce")

    return df


def build_matrix(primary_samples: pd.DataFrame, file_map: dict) -> pd.DataFrame:
    matrices = []

    for _, row in primary_samples.iterrows():
        file_id = row["id"]
        patient_id = row["patient_id"]

        file_path = file_map.get(file_id)

        if file_path is None:
            continue

        df = load_single_file(file_path)

        df = df.set_index("gene_id")
        df.columns = [patient_id]

        matrices.append(df)

    expression_matrix = pd.concat(matrices, axis=1)

    return expression_matrix.T  # patients as rows


def summarize(matrix: pd.DataFrame):
    print("\n Expression Matrix Summary")
    print(f"Shape: {matrix.shape}")
    print(f"Patients: {matrix.shape[0]}")
    print(f"Genes: {matrix.shape[1]}")


def save(matrix: pd.DataFrame):
    out = get_path_from_root("data/interim")
    out.mkdir(parents=True, exist_ok=True)

    file = out / "tcga_blca_expression_matrix_raw_counts.csv"
    matrix.to_csv(file)

    print(f"\n[OK] Saved: {file}")


def main():
    try:
        primary_samples = load_primary_samples()
        file_map = get_file_path_map()

        matrix = build_matrix(primary_samples, file_map)

        summarize(matrix)
        save(matrix)

        return 0

    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())