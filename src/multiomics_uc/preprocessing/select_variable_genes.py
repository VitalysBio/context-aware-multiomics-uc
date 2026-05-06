from __future__ import annotations

import sys
import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_expression() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_expression_log2_cpm_filtered.csv"
    )
    return pd.read_csv(path, index_col=0)


def select_top_variable_genes(
    expression: pd.DataFrame,
    n_genes: int = 5000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_variance = expression.var(axis=0).sort_values(ascending=False)

    selected_genes = gene_variance.head(n_genes).index.tolist()
    selected_expression = expression[selected_genes]

    variance_report = gene_variance.reset_index()
    variance_report.columns = ["gene_id", "variance"]

    return selected_expression, variance_report


def summarize(
    expression: pd.DataFrame,
    selected_expression: pd.DataFrame,
) -> None:
    print("\ Highly Variable Gene Selection Summary ")
    print(f"Input matrix shape: {expression.shape}")
    print(f"Selected matrix shape: {selected_expression.shape}")
    print(f"Selected genes: {selected_expression.shape[1]}")


def save_outputs(
    selected_expression: pd.DataFrame,
    variance_report: pd.DataFrame,
) -> None:
    processed_dir = get_path_from_root("data/processed")
    reports_dir = get_path_from_root("reports/tables")

    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    matrix_file = processed_dir / "tcga_blca_expression_log2_cpm_top5000_hvg.csv"
    report_file = reports_dir / "tcga_blca_gene_variance_report.csv"

    selected_expression.to_csv(matrix_file)
    variance_report.to_csv(report_file, index=False)

    print(f"\n[OK] Saved HVG matrix: {matrix_file}")
    print(f"[OK] Saved variance report: {report_file}")


def main() -> int:
    try:
        expression = load_expression()

        selected_expression, variance_report = select_top_variable_genes(
            expression,
            n_genes=5000,
        )

        summarize(expression, selected_expression)
        save_outputs(selected_expression, variance_report)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())