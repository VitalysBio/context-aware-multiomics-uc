from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_raw_counts() -> pd.DataFrame:
    path = get_path_from_root(
        "data/interim/tcga_blca_expression_matrix_raw_counts.csv"
    )
    return pd.read_csv(path, index_col=0)


def compute_cpm(counts: pd.DataFrame) -> pd.DataFrame:
    library_sizes = counts.sum(axis=1)

    if (library_sizes == 0).any():
        raise ValueError("Some samples have library size equal to zero.")

    cpm = counts.div(library_sizes, axis=0) * 1_000_000
    return cpm


def filter_low_expression_genes(
    cpm: pd.DataFrame,
    min_cpm: float = 1.0,
    min_fraction_samples: float = 0.20,
) -> pd.DataFrame:
    """
    Keep genes with CPM >= min_cpm in at least min_fraction_samples of patients.

    Example:
    min_cpm = 1 and min_fraction_samples = 0.20 means that a gene must have
    CPM >= 1 in at least 20% of patients.
    """
    min_samples = int(np.ceil(cpm.shape[0] * min_fraction_samples))

    keep_genes = (cpm >= min_cpm).sum(axis=0) >= min_samples

    return cpm.loc[:, keep_genes]


def log2_transform(cpm: pd.DataFrame) -> pd.DataFrame:
    return np.log2(cpm + 1)


def summarize(
    raw_counts: pd.DataFrame,
    cpm: pd.DataFrame,
    filtered_cpm: pd.DataFrame,
    log2_cpm: pd.DataFrame,
) -> None:
    print("\n RNA-seq Normalization Summary")
    print(f"Raw matrix shape: {raw_counts.shape}")
    print(f"CPM matrix shape: {cpm.shape}")
    print(f"Filtered CPM matrix shape: {filtered_cpm.shape}")
    print(f"Final log2 CPM matrix shape: {log2_cpm.shape}")

    removed = raw_counts.shape[1] - filtered_cpm.shape[1]
    retained = filtered_cpm.shape[1]

    print(f"\nGenes retained: {retained}")
    print(f"Genes removed: {removed}")


def save_outputs(log2_cpm: pd.DataFrame, filtered_cpm: pd.DataFrame) -> None:
    out_dir = get_path_from_root("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)

    log_file = out_dir / "tcga_blca_expression_log2_cpm_filtered.csv"
    cpm_file = out_dir / "tcga_blca_expression_cpm_filtered.csv"

    log2_cpm.to_csv(log_file)
    filtered_cpm.to_csv(cpm_file)

    print(f"\n[OK] Saved log2 CPM matrix: {log_file}")
    print(f"[OK] Saved filtered CPM matrix: {cpm_file}")


def main() -> int:
    try:
        raw_counts = load_raw_counts()

        cpm = compute_cpm(raw_counts)

        filtered_cpm = filter_low_expression_genes(
            cpm,
            min_cpm=1.0,
            min_fraction_samples=0.20,
        )

        log2_cpm = log2_transform(filtered_cpm)

        summarize(raw_counts, cpm, filtered_cpm, log2_cpm)
        save_outputs(log2_cpm, filtered_cpm)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())