from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import umap

from multiomics_uc.paths import get_path_from_root


def load_expression() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_expression_log2_cpm_top5000_hvg.csv"
    )
    return pd.read_csv(path, index_col=0)


def load_clinical() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_master_clinical_table.csv"
    )
    return pd.read_csv(path)


def scale_expression(expression: pd.DataFrame) -> pd.DataFrame:
    scaler = StandardScaler()
    scaled = scaler.fit_transform(expression)

    return pd.DataFrame(
        scaled,
        index=expression.index,
        columns=expression.columns,
    )


def run_pca(
    expression_scaled: pd.DataFrame,
    n_components: int = 50,
) -> tuple[pd.DataFrame, PCA]:
    pca = PCA(n_components=n_components, random_state=42)
    pcs = pca.fit_transform(expression_scaled)

    pc_cols = [f"PC{i + 1}" for i in range(n_components)]

    pca_df = pd.DataFrame(
        pcs,
        index=expression_scaled.index,
        columns=pc_cols,
    )

    return pca_df, pca


def run_umap(pca_df: pd.DataFrame) -> pd.DataFrame:
    reducer = umap.UMAP(
        n_neighbors=20,
        min_dist=0.20,
        n_components=2,
        metric="euclidean",
        random_state=42,
    )

    embedding = reducer.fit_transform(pca_df)

    umap_df = pd.DataFrame(
        embedding,
        index=pca_df.index,
        columns=["UMAP1", "UMAP2"],
    )

    return umap_df


def merge_with_clinical(
    pca_df: pd.DataFrame,
    umap_df: pd.DataFrame,
    clinical: pd.DataFrame,
) -> pd.DataFrame:
    embeddings = pd.concat([umap_df, pca_df], axis=1).reset_index()
    embeddings = embeddings.rename(columns={"index": "patient_id"})

    keep_cols = [
        "patient_id",
        "sample_barcode",
        "gender",
        "age_at_index",
        "vital_status",
        "survival_time",
        "survival_event",
        "ajcc_pathologic_stage",
        "ajcc_pathologic_t",
        "ajcc_pathologic_n",
        "ajcc_pathologic_m",
    ]

    clinical_subset = clinical[[c for c in keep_cols if c in clinical.columns]]

    merged = embeddings.merge(
        clinical_subset,
        on="patient_id",
        how="left",
        validate="one_to_one",
    )

    return merged


def save_outputs(merged: pd.DataFrame, pca: PCA) -> None:
    processed_dir = get_path_from_root("data/processed")
    reports_dir = get_path_from_root("reports/tables")

    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    embedding_file = processed_dir / "tcga_blca_rnaseq_pca_umap_embeddings.csv"
    variance_file = reports_dir / "tcga_blca_pca_explained_variance.csv"

    merged.to_csv(embedding_file, index=False)

    variance_df = pd.DataFrame(
        {
            "component": [f"PC{i + 1}" for i in range(len(pca.explained_variance_ratio_))],
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_explained_variance": pca.explained_variance_ratio_.cumsum(),
        }
    )
    variance_df.to_csv(variance_file, index=False)

    print(f"\n[OK] Saved embeddings: {embedding_file}")
    print(f"[OK] Saved PCA variance report: {variance_file}")


def plot_umap(merged: pd.DataFrame) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))

    alive = merged["vital_status"] == "Alive"
    dead = merged["vital_status"] == "Dead"

    ax.scatter(
        merged.loc[alive, "UMAP1"],
        merged.loc[alive, "UMAP2"],
        label="Alive",
        alpha=0.75,
        s=35,
    )

    ax.scatter(
        merged.loc[dead, "UMAP1"],
        merged.loc[dead, "UMAP2"],
        label="Dead",
        alpha=0.75,
        s=35,
    )

    ax.set_title("TCGA-BLCA RNA-seq UMAP")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False)

    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_umap_rnaseq_vital_status.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved UMAP figure: {out_file}")


def summarize(merged: pd.DataFrame, pca: PCA) -> None:
    print("\n RNA-seq Embedding Summary")
    print(f"Patients: {merged['patient_id'].nunique()}")
    print(f"Embedding table shape: {merged.shape}")
    print(f"PCA components: {len(pca.explained_variance_ratio_)}")
    print(
        "Cumulative variance explained by 50 PCs: "
        f"{pca.explained_variance_ratio_.sum():.3f}"
    )

    print("\nVital status distribution:")
    print(merged["vital_status"].value_counts(dropna=False))


def main() -> int:
    try:
        expression = load_expression()
        clinical = load_clinical()

        expression_scaled = scale_expression(expression)
        pca_df, pca = run_pca(expression_scaled, n_components=50)
        umap_df = run_umap(pca_df)

        merged = merge_with_clinical(pca_df, umap_df, clinical)

        summarize(merged, pca)
        save_outputs(merged, pca)
        plot_umap(merged)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())