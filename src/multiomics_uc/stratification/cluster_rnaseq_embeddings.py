from __future__ import annotations

import sys
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from multiomics_uc.paths import get_path_from_root


def load_embeddings() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/tcga_blca_rnaseq_pca_umap_embeddings.csv"
    )
    return pd.read_csv(path)


def get_pc_matrix(df: pd.DataFrame) -> pd.DataFrame:
    pc_cols = [col for col in df.columns if col.startswith("PC")]
    if not pc_cols:
        raise ValueError("No PC columns found in embeddings table.")

    return df[pc_cols]


def evaluate_kmeans(pc_matrix: pd.DataFrame, k_values: range) -> pd.DataFrame:
    rows = []

    for k in k_values:
        model = KMeans(
            n_clusters=k,
            random_state=42,
            n_init=50,
        )

        labels = model.fit_predict(pc_matrix)
        sil = silhouette_score(pc_matrix, labels)

        rows.append(
            {
                "k": k,
                "silhouette_score": sil,
                "inertia": model.inertia_,
            }
        )

    return pd.DataFrame(rows)


def assign_clusters(
    df: pd.DataFrame,
    pc_matrix: pd.DataFrame,
    selected_k: int = 4,
) -> pd.DataFrame:
    model = KMeans(
        n_clusters=selected_k,
        random_state=42,
        n_init=50,
    )

    clustered = df.copy()
    clustered["rnaseq_cluster"] = model.fit_predict(pc_matrix)
    clustered["rnaseq_cluster"] = clustered["rnaseq_cluster"].astype(str)

    return clustered


def summarize(clustered: pd.DataFrame, evaluation: pd.DataFrame) -> None:
    print("\n RNA-seq Clustering Summary")
    print("\nKMeans evaluation:")
    print(evaluation.to_string(index=False))

    print("\nCluster sizes:")
    print(clustered["rnaseq_cluster"].value_counts().sort_index())


def save_outputs(clustered: pd.DataFrame, evaluation: pd.DataFrame) -> None:
    processed_dir = get_path_from_root("data/processed")
    reports_dir = get_path_from_root("reports/tables")

    processed_dir.mkdir(parents=True, exist_ok=True)
    reports_dir.mkdir(parents=True, exist_ok=True)

    clustered_file = processed_dir / "tcga_blca_rnaseq_clusters.csv"
    evaluation_file = reports_dir / "tcga_blca_rnaseq_kmeans_evaluation.csv"

    clustered.to_csv(clustered_file, index=False)
    evaluation.to_csv(evaluation_file, index=False)

    print(f"\n[OK] Saved clustered table: {clustered_file}")
    print(f"[OK] Saved clustering evaluation: {evaluation_file}")


def plot_clusters(clustered: pd.DataFrame) -> None:
    figures_dir = get_path_from_root("reports/figures")
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 6))

    for cluster, subset in clustered.groupby("rnaseq_cluster"):
        ax.scatter(
            subset["UMAP1"],
            subset["UMAP2"],
            label=f"Cluster {cluster}",
            alpha=0.75,
            s=35,
        )

    ax.set_title("TCGA-BLCA RNA-seq UMAP by transcriptomic cluster")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(frameon=False, title="Cluster")

    fig.tight_layout()

    out_file = figures_dir / "tcga_blca_umap_rnaseq_clusters.png"
    fig.savefig(out_file, dpi=300)
    plt.close(fig)

    print(f"[OK] Saved cluster UMAP figure: {out_file}")


def main() -> int:
    try:
        embeddings = load_embeddings()
        pc_matrix = get_pc_matrix(embeddings)

        evaluation = evaluate_kmeans(pc_matrix, range(2, 7))

        selected_k = 4
        clustered = assign_clusters(
            embeddings,
            pc_matrix,
            selected_k=selected_k,
        )

        summarize(clustered, evaluation)
        save_outputs(clustered, evaluation)
        plot_clusters(clustered)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())