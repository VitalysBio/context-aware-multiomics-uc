from __future__ import annotations

import argparse
import sys
import pandas as pd

from multiomics_uc.paths import get_path_from_root


PROGRAM_DESCRIPTIONS = {
    "immune_infiltration": "macrophage-rich immune-stromal infiltration",
    "angiogenic_stroma": "angiogenic and vascular stromal remodeling",
    "proliferation": "cell-cycle and proliferative activity",
    "motility_cytoskeleton": "motility and cytoskeletal remodeling",
    "immune_activation": "IFN-associated immune activation",
    "invasion_signaling": "RTK / RHO / MET invasion signaling",
    "differentiation": "epithelial differentiation and ciliary organization",
}


def load_percentiles() -> pd.DataFrame:
    path = get_path_from_root("data/processed/program_percentiles.csv")
    return pd.read_csv(path)


def get_patient_row(df: pd.DataFrame, patient_id: str) -> pd.Series:
    subset = df[df["patient_id"] == patient_id]

    if subset.empty:
        raise ValueError(f"Patient not found: {patient_id}")

    return subset.iloc[0]


def classify_percentile(value: float) -> str:
    if value >= 80:
        return "high"
    if value >= 60:
        return "moderately high"
    if value <= 20:
        return "low"
    if value <= 40:
        return "moderately low"
    return "intermediate"


def build_program_table(row: pd.Series) -> str:
    lines = []

    for program, description in PROGRAM_DESCRIPTIONS.items():
        col = f"{program}_percentile"

        if col not in row.index:
            continue

        value = float(row[col])
        category = classify_percentile(value)

        lines.append(
            f"| {description} | {value:.1f} | {category} |"
        )

    return "\n".join(lines)


def build_narrative(row: pd.Series) -> str:
    high_programs = []
    low_programs = []

    for program, description in PROGRAM_DESCRIPTIONS.items():
        col = f"{program}_percentile"

        if col not in row.index:
            continue

        value = float(row[col])

        if value >= 80:
            high_programs.append(description)

        if value <= 20:
            low_programs.append(description)

    sentences = []

    if high_programs:
        sentences.append(
            "This tumor shows high relative activity in "
            + ", ".join(high_programs)
            + "."
        )

    if low_programs:
        sentences.append(
            "This tumor shows low relative activity in "
            + ", ".join(low_programs)
            + "."
        )

    if not sentences:
        sentences.append(
            "This tumor shows an intermediate biological profile across the major latent programs."
        )

    return " ".join(sentences)


def generate_summary(patient_id: str, row: pd.Series) -> str:
    return f"""# Biological State Summary

## Patient ID

`{patient_id}`

---

## Latent Biological Program Percentiles

| Biological program | Cohort percentile | Category |
|---|---:|---|
{build_program_table(row)}

---

## Automated Biological Interpretation

{build_narrative(row)}

---

## Interpretation Note

Percentiles are computed relative to the TCGA-BLCA cohort used in this project. A high percentile indicates that the patient has stronger expression of that latent biological program compared with most patients in the cohort. A low percentile indicates comparatively reduced activity of that program.
"""


def save_summary(patient_id: str, summary: str) -> None:
    out_dir = get_path_from_root("reports/patient_reports")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / f"{patient_id}_biological_state_summary.md"
    out_file.write_text(summary, encoding="utf-8")

    print(f"[OK] Saved biological state summary: {out_file}")


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--patient-id",
        required=True,
        help="TCGA patient ID",
    )

    args = parser.parse_args()

    try:
        percentiles = load_percentiles()
        row = get_patient_row(percentiles, args.patient_id)

        summary = generate_summary(args.patient_id, row)

        save_summary(args.patient_id, summary)

        print(summary)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())