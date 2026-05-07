from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_annotation() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_gene_annotation_table.csv")
    return pd.read_csv(path)


def marker_files() -> list[Path]:
    marker_dir = get_path_from_root("reports/tables/cluster_markers")
    files = sorted(marker_dir.glob("rnaseq_cluster_*_markers_vs_rest.csv"))

    if not files:
        raise FileNotFoundError(f"No marker files found in: {marker_dir}")

    return files


def annotate_markers(markers: pd.DataFrame, annotation: pd.DataFrame) -> pd.DataFrame:
    annotated = markers.merge(
        annotation,
        on="gene_id",
        how="left",
        validate="many_to_one",
    )

    cols_first = [
        "gene_id",
        "gene_id_clean",
        "gene_name",
        "gene_type",
        "log2_fc",
        "adj_p_value",
        "p_value",
        "mean_cluster",
        "mean_rest",
        "t_statistic",
        "abs_log2_fc",
    ]

    existing_cols = [col for col in cols_first if col in annotated.columns]
    remaining_cols = [col for col in annotated.columns if col not in existing_cols]

    annotated = annotated[existing_cols + remaining_cols]

    return annotated


def summarize(annotated: pd.DataFrame, file: Path) -> None:
    significant = annotated[
        (annotated["adj_p_value"] < 0.05)
        & (annotated["abs_log2_fc"] >= 1)
    ]

    up = significant[significant["log2_fc"] > 0]

    print(f"\n {file.name}")
    print(f"Rows: {len(annotated)}")
    print(f"Annotated gene names: {annotated['gene_name'].notna().sum()}")
    print(f"Significant markers: {len(significant)}")

    print("\nTop upregulated annotated genes:")
    cols = ["gene_name", "gene_type", "log2_fc", "adj_p_value"]
    print(
        up.sort_values(["adj_p_value", "log2_fc"], ascending=[True, False])
        .head(15)[cols]
        .to_string(index=False)
    )


def save_annotated(annotated: pd.DataFrame, source_file: Path) -> None:
    out_dir = get_path_from_root("reports/tables/cluster_markers_annotated")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / source_file.name.replace(
        "_markers_vs_rest.csv",
        "_markers_vs_rest_annotated.csv",
    )

    annotated.to_csv(out_file, index=False)

    print(f"[OK] Saved: {out_file}")


def main() -> int:
    try:
        annotation = load_annotation()
        files = marker_files()

        for file in files:
            markers = pd.read_csv(file)
            annotated = annotate_markers(markers, annotation)
            summarize(annotated, file)
            save_annotated(annotated, file)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())