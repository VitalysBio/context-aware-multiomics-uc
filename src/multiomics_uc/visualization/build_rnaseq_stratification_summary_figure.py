from __future__ import annotations

import sys
from pathlib import Path

from PIL import Image, ImageOps, ImageDraw

from multiomics_uc.paths import get_path_from_root


FIGURES = [
    ("A. RNA-seq UMAP clusters", "tcga_blca_umap_rnaseq_clusters.png"),
    ("B. Overall survival by cluster", "tcga_blca_kaplan_meier_by_rnaseq_cluster.png"),
    ("C. Marker gene heatmap", "tcga_blca_cluster_marker_heatmap.png"),
    ("D. Pathway enrichment summary", "tcga_blca_pathway_enrichment_summary.png"),
]


def load_and_resize(path: Path, width: int) -> Image.Image:
    img = Image.open(path).convert("RGB")

    ratio = width / img.width
    height = int(img.height * ratio)

    return img.resize((width, height), Image.LANCZOS)


def add_panel_title(img: Image.Image, title: str) -> Image.Image:
    title_height = 60
    canvas = Image.new("RGB", (img.width, img.height + title_height), "white")

    draw = ImageDraw.Draw(canvas)
    draw.text((15, 18), title, fill="black")

    canvas.paste(img, (0, title_height))
    return canvas


def build_summary_figure() -> None:
    fig_dir = get_path_from_root("reports/figures")
    output_file = fig_dir / "tcga_blca_rnaseq_stratification_summary.png"

    panel_width = 1000
    panels = []

    for title, filename in FIGURES:
        path = fig_dir / filename

        if not path.exists():
            raise FileNotFoundError(f"Missing figure: {path}")

        img = load_and_resize(path, panel_width)
        img = ImageOps.expand(img, border=20, fill="white")
        img = add_panel_title(img, title)

        panels.append(img)

    row1_height = max(panels[0].height, panels[1].height)
    row2_height = max(panels[2].height, panels[3].height)

    total_width = panel_width * 2 + 80
    total_height = row1_height + row2_height + 80

    canvas = Image.new("RGB", (total_width, total_height), "white")

    canvas.paste(panels[0], (20, 20))
    canvas.paste(panels[1], (panel_width + 60, 20))
    canvas.paste(panels[2], (20, row1_height + 60))
    canvas.paste(panels[3], (panel_width + 60, row1_height + 60))

    canvas.save(output_file, dpi=(300, 300))

    print(f"[OK] Saved summary figure: {output_file}")


def main() -> int:
    try:
        build_summary_figure()
        return 0
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())