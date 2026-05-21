from __future__ import annotations

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from multiomics_uc.paths import get_path_from_root


CANCER_GENES = [
    "FGFR3",
    "ERBB2",
    "EGFR",
    "MYC",
    "CCND1",
    "CDKN2A",
    "RB1",
    "ARID1A",
    "PIK3CA",
    "TP53",
]

EVENT_FILES = {
    "gain_cn_ge_3": "tcga_blca_cnv_gain_cn_ge_3_matrix.csv",
    "amplification_cn_ge_5": "tcga_blca_cnv_amplification_cn_ge_5_matrix.csv",
    "deletion_cn_le_1": "tcga_blca_cnv_deletion_cn_le_1_matrix.csv",
    "deep_deletion_cn_eq_0": "tcga_blca_cnv_deep_deletion_cn_eq_0_matrix.csv",
}


def load_clusters() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_rnaseq_clusters.csv")
    return pd.read_csv(path)[["patient_id", "rnaseq_cluster"]]


def load_annotation() -> pd.DataFrame:
    path = get_path_from_root("reports/tables/cnv/tcga_blca_cnv_gene_annotation.csv")
    return pd.read_csv(path)


def gene_symbol_to_ids(annotation: pd.DataFrame) -> dict[str, str]:
    mapping = {}

    for gene in CANCER_GENES:
        matches = annotation.loc[annotation["gene_name"] == gene, "gene_id"].tolist()
        if matches:
            mapping[gene] = matches[0]

    return mapping


def load_event_matrix(event_file: str) -> pd.DataFrame:
    path = get_path_from_root("data/processed", event_file)
    matrix = pd.read_csv(path, index_col=0)
    matrix.index.name = "patient_id"
    return matrix


def subset_and_rename_genes(
    matrix: pd.DataFrame,
    mapping: dict[str, str],
) -> pd.DataFrame:
    selected = {}

    for symbol, gene_id in mapping.items():
        if gene_id in matrix.columns:
            selected[symbol] = matrix[gene_id]

    out = pd.DataFrame(selected)
    out.index = matrix.index
    out.index.name = "patient_id"

    return out


def merge_clusters(clusters: pd.DataFrame, gene_events: pd.DataFrame) -> pd.DataFrame:
    gene_events = gene_events.reset_index()

    merged = clusters.merge(
        gene_events,
        on="patient_id",
        how="inner",
    )

    merged["rnaseq_cluster"] = merged["rnaseq_cluster"].astype(str)

    return merged


def run_fisher_tests(
    df: pd.DataFrame,
    event_name: str,
    genes: list[str],
) -> pd.DataFrame:
    rows = []

    for cluster in sorted(df["rnaseq_cluster"].unique()):
        in_cluster = df["rnaseq_cluster"] == cluster
        out_cluster = ~in_cluster

        for gene in genes:
            a = ((df[gene] == 1) & in_cluster).sum()
            b = ((df[gene] == 0) & in_cluster).sum()
            c = ((df[gene] == 1) & out_cluster).sum()
            d = ((df[gene] == 0) & out_cluster).sum()

            table = np.array([[a, b], [c, d]])
            odds_ratio, p_value = fisher_exact(table)

            rows.append(
                {
                    "event": event_name,
                    "cluster": cluster,
                    "gene": gene,
                    "freq_cluster": a / (a + b) if (a + b) else np.nan,
                    "freq_other": c / (c + d) if (c + d) else np.nan,
                    "odds_ratio": odds_ratio,
                    "p_value": p_value,
                }
            )

    results = pd.DataFrame(rows)

    results["fdr"] = multipletests(
        results["p_value"],
        method="fdr_bh",
    )[1]

    return results


def build_frequency_table(
    df: pd.DataFrame,
    genes: list[str],
) -> pd.DataFrame:
    return df.groupby("rnaseq_cluster")[genes].mean()


def plot_event_heatmap(
    freq_table: pd.DataFrame,
    event_name: str,
) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(10, 4))

    sns.heatmap(
        freq_table,
        cmap="Reds",
        annot=True,
        fmt=".2f",
        linewidths=0.5,
        vmin=0,
        vmax=1,
    )

    plt.title(f"{event_name} frequency across RNA-seq clusters")
    plt.xlabel("Cancer genes")
    plt.ylabel("RNA-seq cluster")
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)

    plt.tight_layout()

    out_file = figures_dir / f"tcga_blca_cnv_{event_name}_heatmap.png"
    plt.savefig(out_file, dpi=300)
    plt.close()

    print(f"[OK] Saved heatmap: {out_file}")


def save_results(all_results: pd.DataFrame, all_frequencies: dict[str, pd.DataFrame]) -> None:
    out_dir = get_path_from_root("reports/tables/cnv")
    out_dir.mkdir(parents=True, exist_ok=True)

    results_file = out_dir / "cnv_event_enrichment_by_cluster.csv"
    all_results.to_csv(results_file, index=False)

    print(f"[OK] Saved CNV event enrichment results: {results_file}")

    for event_name, freq in all_frequencies.items():
        freq_file = out_dir / f"cnv_{event_name}_frequency_by_cluster.csv"
        freq.to_csv(freq_file)
        print(f"[OK] Saved frequency table: {freq_file}")


def summarize(all_results: pd.DataFrame) -> None:
    print("\nCNV Event Enrichment by RNA-seq Cluster")

    significant = all_results[all_results["fdr"] < 0.05].copy()

    if significant.empty:
        print("\nNo significant CNV event enrichments at FDR < 0.05")
    else:
        significant = significant.sort_values("fdr")
        print("\nTop significant CNV event enrichments:")
        print(
            significant[
                [
                    "event",
                    "cluster",
                    "gene",
                    "freq_cluster",
                    "freq_other",
                    "odds_ratio",
                    "fdr",
                ]
            ]
            .head(30)
            .to_string(index=False)
        )


def main() -> int:
    try:
        clusters = load_clusters()
        annotation = load_annotation()
        mapping = gene_symbol_to_ids(annotation)

        print("\nCancer genes mapped:")
        print(mapping)

        all_results = []
        all_frequencies = {}

        for event_name, event_file in EVENT_FILES.items():
            matrix = load_event_matrix(event_file)
            gene_events = subset_and_rename_genes(matrix, mapping)

            df = merge_clusters(clusters, gene_events)
            genes = list(gene_events.columns)

            results = run_fisher_tests(df, event_name, genes)
            freq_table = build_frequency_table(df, genes)

            all_results.append(results)
            all_frequencies[event_name] = freq_table

            plot_event_heatmap(freq_table, event_name)

        combined = pd.concat(all_results, axis=0, ignore_index=True)

        summarize(combined)
        save_results(combined, all_frequencies)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())