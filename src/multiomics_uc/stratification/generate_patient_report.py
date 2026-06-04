from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from multiomics_uc.paths import get_path_from_root


CLUSTER_PROFILES = {
    "0": {
        "name": "Stromal / ECM-rich",
        "biology": [
            "Extracellular matrix remodeling",
            "Collagen organization",
            "EMT-associated signaling",
            "Stromal-rich tumor microenvironment",
        ],
        "interpretation": (
            "This profile suggests a tumor state dominated by stromal remodeling, "
            "extracellular matrix organization, and EMT-like programs."
        ),
    },
    "1": {
        "name": "Luminal-like metabolic",
        "biology": [
            "Fatty acid metabolism",
            "Lipid metabolism",
            "Estrogen response",
            "Epithelial differentiation",
        ],
        "interpretation": (
            "This profile suggests a more differentiated luminal-like tumor state "
            "with metabolic and epithelial organization programs."
        ),
    },
    "2": {
        "name": "Basal / squamous inflammatory",
        "biology": [
            "Keratinization",
            "Epidermis development",
            "TNF-alpha / NF-kB signaling",
            "Inflammatory response",
            "High CDKN2A deletion tendency",
        ],
        "interpretation": (
            "This profile suggests a basal/squamous tumor state characterized by "
            "keratinization, inflammatory signaling, and frequent CDKN2A loss."
        ),
    },
    "3": {
        "name": "Immune-microenvironment rich",
        "biology": [
            "Macrophage-rich immune infiltration",
            "Complement activation",
            "Innate immune response",
            "IFN-associated signaling",
            "EMT / microenvironment remodeling",
        ],
        "interpretation": (
            "This profile suggests a tumor state enriched in immune and "
            "microenvironmental signals, including myeloid/macrophage-associated "
            "programs and stromal remodeling."
        ),
    },
}


PROGRAM_DESCRIPTIONS = {
    "immune_infiltration": "macrophage-rich immune-stromal infiltration",
    "angiogenic_stroma": "angiogenic and vascular stromal remodeling",
    "proliferation": "cell-cycle and proliferative activity",
    "motility_cytoskeleton": "motility and cytoskeletal remodeling",
    "immune_activation": "IFN-associated immune activation",
    "invasion_signaling": "RTK / RHO / MET invasion signaling",
    "differentiation": "epithelial differentiation and ciliary organization",
}


DRIVER_MUTATION_COLUMNS = [
    "TP53_mut",
    "FGFR3_mut",
    "RB1_mut",
    "KDM6A_mut",
    "ARID1A_mut",
    "PIK3CA_mut",
]


def load_predictions() -> pd.DataFrame:
    path = get_path_from_root(
        "reports/tables/stratification/patient_stratification_cv_predictions.csv"
    )
    return pd.read_csv(path)


def load_features() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/integrated_patient_feature_table_clean.csv"
    )
    return pd.read_csv(path)


def load_percentiles() -> pd.DataFrame:
    path = get_path_from_root("data/processed/program_percentiles.csv")
    return pd.read_csv(path)


def get_patient_row(df: pd.DataFrame, patient_id: str) -> pd.Series:
    subset = df[df["patient_id"] == patient_id]

    if subset.empty:
        raise ValueError(f"Patient not found: {patient_id}")

    return subset.iloc[0]


def classify_uncertainty(row: pd.Series) -> str:
    confidence = float(row["confidence_score"])

    if confidence >= 0.55:
        return "High confidence"
    if confidence >= 0.30:
        return "Moderate confidence"
    return "Low confidence / ambiguous assignment"


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

def interpret_mutation_burden_percentile(value: float) -> str:
    if pd.isna(value):
        return "Mutation burden percentile is not available."

    if value >= 80:
        return (
            "This tumor has a high mutation burden relative to the TCGA-BLCA cohort, "
            "which may reflect increased genomic instability and potentially higher neoantigen load."
        )

    if value >= 60:
        return (
            "This tumor has a moderately elevated mutation burden relative to the TCGA-BLCA cohort."
        )

    if value <= 20:
        return (
            "This tumor has a low mutation burden relative to the TCGA-BLCA cohort."
        )

    if value <= 40:
        return (
            "This tumor has a moderately low mutation burden relative to the TCGA-BLCA cohort."
        )

    return (
        "This tumor has an intermediate mutation burden relative to the TCGA-BLCA cohort."
    )

def get_probability_table(row: pd.Series) -> str:
    prob_cols = [
        col for col in row.index
        if col.startswith("prob_cluster_")
    ]

    lines = []

    for col in sorted(prob_cols):
        cluster = col.replace("prob_cluster_", "")
        name = CLUSTER_PROFILES.get(cluster, {}).get("name", "Unknown")
        probability = float(row[col])

        lines.append(
            f"| Cluster {cluster} | {name} | {probability:.3f} |"
        )

    return "\n".join(lines)


def get_mutation_summary(feature_row: pd.Series) -> str:
    lines = []

    for col in DRIVER_MUTATION_COLUMNS:
        if col not in feature_row.index:
            continue

        gene = col.replace("_mut", "")
        value = feature_row[col]

        if pd.isna(value):
            status = "Not available"
        elif int(value) == 1:
            status = "Mutated"
        else:
            status = "Wild-type / not mutated"

        lines.append(f"| {gene} | {status} |")

    return "\n".join(lines)


def build_program_percentile_table(percentile_row: pd.Series) -> str:
    lines = []

    for program, description in PROGRAM_DESCRIPTIONS.items():
        col = f"{program}_percentile"

        if col not in percentile_row.index:
            continue

        value = float(percentile_row[col])
        category = classify_percentile(value)

        lines.append(
            f"| {description} | {value:.1f} | {category} |"
        )

    return "\n".join(lines)


def build_biological_narrative(percentile_row: pd.Series) -> str:
    high_programs = []
    low_programs = []
    moderate_high_programs = []

    for program, description in PROGRAM_DESCRIPTIONS.items():
        col = f"{program}_percentile"

        if col not in percentile_row.index:
            continue

        value = float(percentile_row[col])

        if value >= 80:
            high_programs.append(description)
        elif value >= 60:
            moderate_high_programs.append(description)

        if value <= 20:
            low_programs.append(description)

    sentences = []

    if high_programs:
        sentences.append(
            "The strongest biological signals in this patient are "
            + ", ".join(high_programs)
            + "."
        )

    if moderate_high_programs:
        sentences.append(
            "Moderately elevated programs include "
            + ", ".join(moderate_high_programs)
            + "."
        )

    if low_programs:
        sentences.append(
            "Low relative programs include "
            + ", ".join(low_programs)
            + "."
        )

    if not sentences:
        sentences.append(
            "This patient shows an intermediate profile across the major latent biological programs."
        )

    return " ".join(sentences)


def build_similarity_comment(prediction_row: pd.Series) -> str:
    prob_cols = [
        col for col in prediction_row.index
        if col.startswith("prob_cluster_")
    ]

    probs = {
        col.replace("prob_cluster_", ""): float(prediction_row[col])
        for col in prob_cols
    }

    sorted_probs = sorted(
        probs.items(),
        key=lambda x: x[1],
        reverse=True,
    )

    top_cluster, top_prob = sorted_probs[0]
    second_cluster, second_prob = sorted_probs[1]

    top_name = CLUSTER_PROFILES[top_cluster]["name"]
    second_name = CLUSTER_PROFILES[second_cluster]["name"]

    if second_prob >= 0.20:
        return (
            f"The tumor is primarily assigned to Cluster {top_cluster} "
            f"({top_name}), but it also shows partial similarity to Cluster "
            f"{second_cluster} ({second_name})."
        )

    return (
        f"The tumor is predominantly assigned to Cluster {top_cluster} "
        f"({top_name}), with limited similarity to other clusters."
    )


def generate_report(
    patient_id: str,
    prediction_row: pd.Series,
    feature_row: pd.Series,
    percentile_row: pd.Series,
) -> str:
    predicted_cluster = str(prediction_row["predicted_cluster"])
    true_cluster = str(prediction_row["rnaseq_cluster"])

    profile = CLUSTER_PROFILES[predicted_cluster]
    uncertainty_label = classify_uncertainty(prediction_row)

    mutation_burden = feature_row.get("mutation_burden", "NA")

    mutation_burden_percentile = percentile_row.get(
        "mutation_burden_percentile",
        None,
    )

    mutation_burden_interpretation = interpret_mutation_burden_percentile(
        mutation_burden_percentile
    )

    report = f"""# Patient Stratification Report

## Patient ID

`{patient_id}`

---

## Predicted Tumor State

**Predicted cluster:** Cluster {predicted_cluster}  
**Predicted biological label:** {profile["name"]}  
**Cross-validated reference cluster:** Cluster {true_cluster}  
**Assignment confidence:** {uncertainty_label}  

| Metric | Value |
|---|---:|
| Max probability | {float(prediction_row["max_probability"]):.3f} |
| Uncertainty entropy | {float(prediction_row["uncertainty_entropy"]):.3f} |
| Confidence score | {float(prediction_row["confidence_score"]):.3f} |

---

## Cluster Probability Profile

| Cluster | Biological label | Probability |
|---|---|---:|
{get_probability_table(prediction_row)}

{build_similarity_comment(prediction_row)}

---

## Cluster-Level Biological Interpretation

{profile["interpretation"]}

Main biological programs associated with this predicted cluster:

{chr(10).join([f"- {item}" for item in profile["biology"]])}

---

## Patient-Specific Latent Biological Program Percentiles

| Biological program | Cohort percentile | Category |
|---|---:|---|
{build_program_percentile_table(percentile_row)}

---

## Patient-Specific Biological State Interpretation

{build_biological_narrative(percentile_row)}

Percentiles are computed relative to the TCGA-BLCA cohort used in this project. High percentiles indicate that the patient has stronger activity of that latent biological program compared with most patients in the cohort.

---

## Driver Mutation Profile

| Gene | Status |
|---|---|
{get_mutation_summary(feature_row)}

---

## Mutation Burden Context

| Feature | Value |
|---|---:|
| Mutation burden | {mutation_burden} |
| Mutation burden percentile | {float(mutation_burden_percentile):.1f} |

{mutation_burden_interpretation}

---

## Translational Interpretation

This report summarizes a context-aware biological assignment derived from transcriptomic latent programs, driver mutation features, and clinical information. The predicted tumor state should be interpreted as a probabilistic biological profile rather than a deterministic diagnostic label.

A high-confidence assignment suggests that the patient strongly resembles one learned TCGA-BLCA transcriptomic state. Lower-confidence assignments may indicate biological ambiguity, mixed tumor programs, or transitional states between clusters.

---

## Disclaimer

This report is generated for research and portfolio demonstration purposes only. It is not intended for clinical diagnosis, treatment selection, or medical decision-making.
"""

    return report


def save_report(patient_id: str, report: str) -> Path:
    out_dir = get_path_from_root("reports/patient_reports")
    out_dir.mkdir(parents=True, exist_ok=True)

    safe_id = patient_id.replace("/", "_")
    out_file = out_dir / f"{safe_id}_patient_stratification_report_v2.md"

    out_file.write_text(report, encoding="utf-8")

    return out_file


def main() -> int:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--patient-id",
        required=True,
        help="TCGA patient ID, e.g. TCGA-2F-A9KO",
    )

    args = parser.parse_args()

    try:
        predictions = load_predictions()
        features = load_features()
        percentiles = load_percentiles()

        prediction_row = get_patient_row(
            predictions,
            args.patient_id,
        )

        feature_row = get_patient_row(
            features,
            args.patient_id,
        )

        percentile_row = get_patient_row(
            percentiles,
            args.patient_id,
        )

        report = generate_report(
            args.patient_id,
            prediction_row,
            feature_row,
            percentile_row,
        )

        out_file = save_report(
            args.patient_id,
            report,
        )

        print(f"[OK] Saved patient report v2: {out_file}")

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())