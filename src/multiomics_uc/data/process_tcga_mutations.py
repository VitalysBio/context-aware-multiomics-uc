from __future__ import annotations

import sys
import gzip
from pathlib import Path

import pandas as pd

from multiomics_uc.paths import get_path_from_root


NON_SILENT_VARIANTS = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Nonstop_Mutation",
    "Translation_Start_Site",
    "In_Frame_Del",
    "In_Frame_Ins",
}


MIN_MUTATED_PATIENTS = 10


def get_maf_files() -> list[Path]:
    maf_dir = get_path_from_root(
        "data/raw/tcga_blca/mutations/maf"
    )

    maf_files = sorted(maf_dir.rglob("*.maf.gz"))

    if not maf_files:
        raise FileNotFoundError("No MAF files found.")

    return maf_files


def load_single_maf(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        df = pd.read_csv(
            handle,
            sep="\t",
            comment="#",
            low_memory=False,
        )

    return df


def extract_patient_id(barcode: str) -> str:
    if pd.isna(barcode):
        return None

    parts = str(barcode).split("-")

    if len(parts) < 3:
        return None

    return "-".join(parts[:3])


def process_maf(df: pd.DataFrame) -> pd.DataFrame:
    required_cols = [
        "Hugo_Symbol",
        "Variant_Classification",
        "Tumor_Sample_Barcode",
        "Variant_Type",
    ]

    missing = [c for c in required_cols if c not in df.columns]

    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    df = df[required_cols].copy()

    df = df[
        df["Variant_Classification"].isin(NON_SILENT_VARIANTS)
    ].copy()

    df["patient_id"] = df["Tumor_Sample_Barcode"].apply(
        extract_patient_id
    )

    df = df.dropna(subset=["patient_id", "Hugo_Symbol"])

    return df


def concatenate_mutations(maf_files: list[Path]) -> pd.DataFrame:
    all_dfs = []

    for idx, maf_file in enumerate(maf_files, start=1):
        try:
            df = load_single_maf(maf_file)
            df = process_maf(df)

            all_dfs.append(df)

            if idx % 50 == 0:
                print(f"Processed {idx}/{len(maf_files)} MAF files")

        except Exception as exc:
            print(f"[WARNING] Failed: {maf_file.name} -> {exc}")

    combined = pd.concat(all_dfs, axis=0, ignore_index=True)

    return combined


def compute_mutation_burden(df: pd.DataFrame) -> pd.DataFrame:
    burden = (
        df.groupby("patient_id")
        .size()
        .reset_index(name="mutation_burden")
    )

    return burden


def select_recurrent_genes(df: pd.DataFrame) -> pd.DataFrame:
    recurrent = (
        df.groupby("Hugo_Symbol")["patient_id"]
        .nunique()
        .reset_index(name="mutated_patients")
        .sort_values(
            "mutated_patients",
            ascending=False,
        )
    )

    recurrent["mutation_frequency"] = (
        recurrent["mutated_patients"]
        / df["patient_id"].nunique()
    )

    return recurrent


def build_binary_matrix(
    df: pd.DataFrame,
    recurrent_genes: pd.DataFrame,
) -> pd.DataFrame:
    selected_genes = recurrent_genes[
        recurrent_genes["mutated_patients"]
        >= MIN_MUTATED_PATIENTS
    ]["Hugo_Symbol"].tolist()

    filtered = df[
        df["Hugo_Symbol"].isin(selected_genes)
    ].copy()

    filtered["mutated"] = 1

    matrix = (
        filtered.pivot_table(
            index="patient_id",
            columns="Hugo_Symbol",
            values="mutated",
            aggfunc="max",
            fill_value=0,
        )
    )

    matrix = matrix.astype(int)

    return matrix


def save_outputs(
    mutation_df: pd.DataFrame,
    burden_df: pd.DataFrame,
    recurrent_df: pd.DataFrame,
    binary_matrix: pd.DataFrame,
) -> None:
    processed_dir = get_path_from_root("data/processed")
    reports_dir = get_path_from_root("reports/tables/mutations")

    reports_dir.mkdir(parents=True, exist_ok=True)

    mutation_df.to_csv(
        processed_dir / "tcga_blca_mutation_table_long.csv",
        index=False,
    )

    burden_df.to_csv(
        processed_dir / "tcga_blca_mutation_burden.csv",
        index=False,
    )

    binary_matrix.to_csv(
        processed_dir / "tcga_blca_mutation_binary_matrix.csv"
    )

    recurrent_df.to_csv(
        reports_dir / "recurrent_driver_genes.csv",
        index=False,
    )

    print("\n[OK] Mutation outputs saved")


def summarize(
    mutation_df: pd.DataFrame,
    burden_df: pd.DataFrame,
    recurrent_df: pd.DataFrame,
    binary_matrix: pd.DataFrame,
) -> None:
    print("\n Mutation Processing Summary")

    print(f"Patients with mutations: {mutation_df['patient_id'].nunique()}")
    print(f"Total non-silent mutations: {len(mutation_df)}")

    print("\nTop recurrent genes:")
    print(recurrent_df.head(15).to_string(index=False))

    print("\nMutation burden summary:")
    print(burden_df["mutation_burden"].describe())

    print("\nBinary matrix shape:")
    print(binary_matrix.shape)


def main() -> int:
    try:
        maf_files = get_maf_files()

        print(f"Found {len(maf_files)} MAF files")

        mutation_df = concatenate_mutations(maf_files)

        burden_df = compute_mutation_burden(mutation_df)

        recurrent_df = select_recurrent_genes(mutation_df)

        binary_matrix = build_binary_matrix(
            mutation_df,
            recurrent_df,
        )

        summarize(
            mutation_df,
            burden_df,
            recurrent_df,
            binary_matrix,
        )

        save_outputs(
            mutation_df,
            burden_df,
            recurrent_df,
            binary_matrix,
        )

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())