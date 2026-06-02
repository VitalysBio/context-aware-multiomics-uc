from __future__ import annotations

import argparse
import sys
import time

import pandas as pd
import gseapy as gp

from multiomics_uc.paths import get_path_from_root


GENE_SET_LIBRARIES = [
    "MSigDB_Hallmark_2020",
    "Reactome_2022",
    "GO_Biological_Process_2023",
]


def load_gene_list(pc: str, direction: str) -> list[str]:
    path = get_path_from_root(
        "reports/tables/pca_interpretation",
        f"{pc}_top_{direction}_genes.csv",
    )

    df = pd.read_csv(path)

    df = df[
        (df["gene_name"].notna())
        & (df["gene_type"] == "protein_coding")
        & (df["fdr"] < 0.05)
    ].copy()

    genes = df["gene_name"].drop_duplicates().astype(str).tolist()

    return genes


def run_enrichment(genes: list[str]) -> pd.DataFrame:
    results = []

    for library in GENE_SET_LIBRARIES:
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

        time.sleep(5)

    return pd.concat(results, ignore_index=True)


def save_results(enrichment: pd.DataFrame, pc: str, direction: str) -> None:
    out_dir = get_path_from_root("reports/tables/pca_interpretation")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"{pc}_{direction}_pathway_enrichment.csv"
    enrichment.to_csv(out_file, index=False)

    print(f"[OK] Saved: {out_file}")


def summarize(enrichment: pd.DataFrame, pc: str, direction: str, n_genes: int) -> None:
    print(f"\n {pc} {direction} gene enrichment")
    print(f"Input genes: {n_genes}")

    top = enrichment.sort_values("Adjusted P-value").head(15)

    cols = [
        "gene_set_library",
        "Term",
        "Adjusted P-value",
        "Odds Ratio",
        "Combined Score",
    ]

    cols = [c for c in cols if c in top.columns]

    print(top[cols].to_string(index=False))


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--pc",
        type=str,
        required=True,
        help="PCA component to analyze, e.g. PC8 or PC9.",
    )

    parser.add_argument(
        "--direction",
        type=str,
        required=True,
        choices=["positive", "negative"],
        help="Gene direction to analyze.",
    )

    args = parser.parse_args()

    try:
        genes = load_gene_list(args.pc, args.direction)

        if len(genes) < 5:
            raise ValueError(f"Too few genes for {args.pc} {args.direction}")

        enrichment = run_enrichment(genes)

        summarize(enrichment, args.pc, args.direction, len(genes))
        save_results(enrichment, args.pc, args.direction)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())