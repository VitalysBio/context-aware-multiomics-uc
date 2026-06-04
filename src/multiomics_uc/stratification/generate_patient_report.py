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
            "extracellular matrix organization, and EMT-like programs. This state was "
            "associated with poorer prognosis in the transcriptomic survival analysis."
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
            "This profile suggests a more differentiated luminal-like tumor state with "
            "metabolic and epithelial organization programs. This cluster showed a more "
            "favorable survival profile compared with stromal and basal-like groups."
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
            "This profile suggests a basal/squamous tumor state characterized by keratinization, "
            "inflammatory signaling, and frequent CDKN2A loss. This group was associated with "
            "poorer prognosis and a more aggressive transcriptomic phenotype."
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
            "This profile suggests a tumor state enriched in immune and microenvironmental "
            "signals, including myeloid/macrophage-associated programs and stromal remodeling. "
            "This cluster may represent a biologically mixed immune-stromal phenotype."
        ),
    },
}


LATENT_PROGRAMS = {
    "PC1": "Macrophage-rich immune-stromal infiltration",
    "PC2": "Angiogenic / vascular stromal microenvironment versus basal-squamous axis",
    "PC3": "Cell-cycle, mitosis and proliferation",
    "PC4": "Motility and cytoskeletal remodeling",
    "PC6": "IFN-driven immune activation versus EMT / ECM remodeling",
    "PC8": "RTK / RHO / MET invasion signaling",
    "PC9": "Differentiation, cilium organization and epithelial structure",
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


def get_probability_table(row: pd.Series) -> str:
    prob_cols = [
        col for col in row.index
        if col.startswith("prob_cluster_")
    ]

    lines = []

    for col in sorted(prob_cols):
        cluster = col.replace("prob_cluster_", "")
        profile = CLUSTER_PROFILES.get(cluster, {})
        name = profile.get("name", "Unknown")
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


def get_latent_program_summary(feature_row: pd.Series) -> str:
    pcs = ["PC1", "PC2", "PC3", "PC4", "PC6", "PC8", "PC9"]

    lines = []

    for pc in pcs:
        if pc not in feature_row.index:
            continue

        value = float(feature_row[pc])
        meaning = LATENT_PROGRAMS[pc]

        direction = "High" if value > 0 else "Low"

        lines.append(
            f"| {pc} | {value:.3f} | {direction} | {meaning} |"
        )

    return "\n".join(lines)


def generate_report(
    patient_id: str,
    prediction_row: pd.Series,
    feature_row: pd.Series,
) -> str:
    predicted_cluster = str(prediction_row["predicted_cluster"])
    true_cluster = str(prediction_row["rnaseq_cluster"])

    profile = CLUSTER_PROFILES[predicted_cluster]

    uncertainty_label = classify_uncertainty(prediction_row)

    mutation_burden = feature_row.get("mutation_burden", "NA")

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

---

## Biological Interpretation

{profile["interpretation"]}

Main biological programs associated with this predicted state:

{chr(10).join([f"- {item}" for item in profile["biology"]])}

---

## Latent Transcriptomic Program Profile

| Program | Patient score | Direction | Biological interpretation |
|---|---:|---|---|
{get_latent_program_summary(feature_row)}

---

## Driver Mutation Profile

| Gene | Status |
|---|---|
{get_mutation_summary(feature_row)}

---

## Mutation Burden

| Feature | Value |
|---|---:|
| Mutation burden | {mutation_burden} |

---

## Translational Interpretation

This report summarizes a context-aware biological assignment derived from transcriptomic latent programs, driver mutation features, and clinical information. The predicted tumor state should be interpreted as a probabilistic biological profile rather than a deterministic diagnostic label.

High-confidence assignments suggest that the patient strongly resembles one of the learned TCGA-BLCA transcriptomic states. Low-confidence assignments may indicate biological ambiguity, mixed tumor programs, or transitional states between clusters.

---

## Disclaimer

This report is generated for research and portfolio demonstration purposes only. It is not intended for clinical diagnosis, treatment selection, or medical decision-making.
"""

    return report


def save_report(patient_id: str, report: str) -> Path:
    out_dir = get_path_from_root("reports/patient_reports")
    out_dir.mkdir(parents=True, exist_ok=True)

    safe_id = patient_id.replace("/", "_")
    out_file = out_dir / f"{safe_id}_patient_stratification_report.md"

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

        prediction_row = get_patient_row(
            predictions,
            args.patient_id,
        )

        feature_row = get_patient_row(
            features,
            args.patient_id,
        )

        report = generate_report(
            args.patient_id,
            prediction_row,
            feature_row,
        )

        out_file = save_report(
            args.patient_id,
            report,
        )

        print(f"[OK] Saved patient report: {out_file}")

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())