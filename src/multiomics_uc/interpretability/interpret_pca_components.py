from __future__ import annotations

import sys
import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

from multiomics_uc.paths import get_path_from_root


PCS_TO_INTERPRET = ["PC8", "PC9"]
TOP_N = 300


def load_expression() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_expression_log2_cpm_filtered.csv"
    )
    return pd.read_csv(path, index_col=0)


def load_embeddings() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_rnaseq_pca_umap_embeddings.csv"
    )
    return pd.read_csv(path).set_index("patient_id")


def load_annotation() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_gene_annotation_table.csv")
    return pd.read_csv(path)


def align_data(
    expression: pd.DataFrame,
    embeddings: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    common = expression.index.intersection(embeddings.index)

    if len(common) == 0:
        raise ValueError("No common patients between expression and embeddings.")

    return expression.loc[common], embeddings.loc[common]


def correlate_genes_with_pc(
    expression: pd.DataFrame,
    pc_values: pd.Series,
) -> pd.DataFrame:
    rows = []

    for gene in expression.columns:
        r, p = pearsonr(expression[gene], pc_values)

        rows.append(
            {
                "gene_id": gene,
                "correlation": r,
                "p_value": p,
                "abs_correlation": abs(r),
            }
        )

    corr_df = pd.DataFrame(rows)

    corr_df["fdr"] = multipletests(
        corr_df["p_value"],
        method="fdr_bh",
    )[1]

    corr_df = corr_df.sort_values(
        "abs_correlation",
        ascending=False,
    )

    return corr_df


def annotate_correlations(
    corr: pd.DataFrame,
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    return corr.merge(
        annotation,
        on="gene_id",
        how="left",
    )


def save_pc_results(pc: str, annotated: pd.DataFrame) -> None:
    out_dir = get_path_from_root("reports/tables/pca_interpretation")
    out_dir.mkdir(parents=True, exist_ok=True)

    full_file = out_dir / f"{pc}_gene_correlations.csv"
    pos_file = out_dir / f"{pc}_top_positive_genes.csv"
    neg_file = out_dir / f"{pc}_top_negative_genes.csv"

    annotated.to_csv(full_file, index=False)

    annotated.sort_values("correlation", ascending=False).head(TOP_N).to_csv(
        pos_file,
        index=False,
    )

    annotated.sort_values("correlation", ascending=True).head(TOP_N).to_csv(
        neg_file,
        index=False,
    )

    print(f"[OK] Saved PC interpretation files for {pc}")


def summarize(pc: str, annotated: pd.DataFrame) -> None:
    print(f"\n {pc} interpretation ")

    significant = annotated[
        (annotated["fdr"] < 0.05)
        & (annotated["abs_correlation"] >= 0.30)
    ]

    print(f"Total genes tested: {len(annotated)}")
    print(
        "Significant genes with FDR < 0.05 and |correlation| >= 0.30: "
        f"{len(significant)}"
    )

    cols = [
        "gene_name",
        "gene_type",
        "correlation",
        "p_value",
        "fdr",
    ]

    print("\nTop positively correlated genes:")
    print(
        annotated.sort_values("correlation", ascending=False)
        .head(15)[cols]
        .to_string(index=False)
    )

    print("\nTop negatively correlated genes:")
    print(
        annotated.sort_values("correlation", ascending=True)
        .head(15)[cols]
        .to_string(index=False)
    )


def main() -> int:
    try:
        expression = load_expression()
        embeddings = load_embeddings()
        annotation = load_annotation()

        expression, embeddings = align_data(expression, embeddings)

        print("\n PCA Component Biological Interpretation")
        print(f"Patients: {expression.shape[0]}")
        print(f"Genes tested: {expression.shape[1]}")
        print(f"PCs interpreted: {PCS_TO_INTERPRET}")

        for pc in PCS_TO_INTERPRET:
            if pc not in embeddings.columns:
                raise ValueError(f"{pc} not found in embeddings table.")

            corr = correlate_genes_with_pc(
                expression,
                embeddings[pc],
            )

            annotated = annotate_correlations(
                corr,
                annotation,
            )

            summarize(pc, annotated)
            save_pc_results(pc, annotated)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())