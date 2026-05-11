from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import gseapy as gp
import argparse
import time

from multiomics_uc.paths import get_path_from_root


GENE_SET_LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
    "GO_Biological_Process_2023",
]


def marker_files() -> list[Path]:
    marker_dir = get_path_from_root("reports/tables/cluster_markers_annotated")
    files = sorted(marker_dir.glob("rnaseq_cluster_*_markers_vs_rest_annotated.csv"))

    if not files:
        raise FileNotFoundError(f"No annotated marker files found in: {marker_dir}")

    return files


def extract_cluster_label(file_path: Path) -> str:
    # Example:
    # rnaseq_cluster_0_markers_vs_rest_annotated.csv
    return file_path.name.split("_")[2]


def load_upregulated_genes(
    file_path: Path,
    fdr_threshold: float = 0.05,
    log2fc_threshold: float = 1.0,
    protein_coding_only: bool = True,
) -> list[str]:
    df = pd.read_csv(file_path)

    df = df[
        (df["adj_p_value"] < fdr_threshold)
        & (df["log2_fc"] >= log2fc_threshold)
        & df["gene_name"].notna()
    ].copy()

    if protein_coding_only:
        df = df[df["gene_type"] == "protein_coding"]

    genes = (
        df["gene_name"]
        .astype(str)
        .drop_duplicates()
        .tolist()
    )

    return genes


def run_enrichment(
    genes: list[str],
    gene_sets: list[str],
) -> pd.DataFrame:
    if len(genes) < 5:
        raise ValueError("Too few genes for enrichment analysis.")

    results = []

    for library in gene_sets:
        enr = gp.enrichr(
            gene_list=genes,
            gene_sets=library,
            organism="human",
            outdir=None,
            cutoff=0.5,
            no_plot=True,
        )

        res = enr.results.copy()
        res["gene_set_library"] = library
        results.append(res)

    combined = pd.concat(results, axis=0, ignore_index=True)

    return combined


def save_enrichment(
    enrichment: pd.DataFrame,
    cluster_label: str,
) -> None:
    out_dir = get_path_from_root("reports/tables/pathway_enrichment")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"rnaseq_cluster_{cluster_label}_pathway_enrichment.csv"
    enrichment.to_csv(out_file, index=False)

    print(f"[OK] Saved enrichment for cluster {cluster_label}: {out_file}")


def summarize(
    enrichment: pd.DataFrame,
    cluster_label: str,
    n_genes: int,
) -> None:
    print(f"\n=== Cluster {cluster_label} pathway enrichment ===")
    print(f"Input upregulated protein-coding genes: {n_genes}")

    if "Adjusted P-value" not in enrichment.columns:
        print("[WARN] Adjusted P-value column not found.")
        print(enrichment.head().to_string(index=False))
        return

    top = enrichment.sort_values("Adjusted P-value").head(15)

    cols = [
        "gene_set_library",
        "Term",
        "Adjusted P-value",
        "Odds Ratio",
        "Combined Score",
    ]
    cols = [c for c in cols if c in top.columns]

    print("\nTop enriched pathways:")
    print(top[cols].to_string(index=False))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cluster",
        type=str,
        default=None,
        help="Run enrichment for one cluster only, e.g. 0, 1, 2, or 3.",
    )
    args = parser.parse_args()

    try:
        files = marker_files()

        if args.cluster is not None:
            files = [
                f for f in files
                if extract_cluster_label(f) == str(args.cluster)
            ]

            if not files:
                raise ValueError(f"No marker file found for cluster {args.cluster}")

        print("\n Pathway Enrichment by RNA-seq Cluster")

        for file_path in files:
            cluster_label = extract_cluster_label(file_path)

            genes = load_upregulated_genes(
                file_path,
                fdr_threshold=0.05,
                log2fc_threshold=1.0,
                protein_coding_only=True,
            )

            enrichment = run_enrichment(
                genes,
                gene_sets=GENE_SET_LIBRARIES,
            )

            summarize(enrichment, cluster_label, len(genes))
            save_enrichment(enrichment, cluster_label)

            time.sleep(10)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())