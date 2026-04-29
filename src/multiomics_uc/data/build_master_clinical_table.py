from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_primary_samples() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_primary_samples.csv")
    return pd.read_csv(path)


def load_clinical_json() -> list[dict]:
    path = get_path_from_root("data/external/tcga_blca/metadata/clinical.json")

    if not path.exists():
        raise FileNotFoundError(f"Clinical JSON not found: {path}")

    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def safe_get(dictionary: dict, key: str, default=None):
    if not isinstance(dictionary, dict):
        return default
    return dictionary.get(key, default)


def first_item(items: list):
    if isinstance(items, list) and len(items) > 0:
        return items[0]
    return {}


def parse_clinical_record(record: dict) -> dict:
    demographic = safe_get(record, "demographic", {})
    diagnoses = safe_get(record, "diagnoses", [])
    diagnosis = first_item(diagnoses)

    patient_id = record.get("submitter_id")

    vital_status = safe_get(demographic, "vital_status")
    days_to_death = safe_get(demographic, "days_to_death")
    days_to_last_follow_up = safe_get(diagnosis, "days_to_last_follow_up")

    survival_time = days_to_death if pd.notnull(days_to_death) else days_to_last_follow_up
    survival_event = 1 if vital_status == "Dead" else 0

    return {
        "patient_id": patient_id,
        "case_id": record.get("case_id"),
        "disease_type": record.get("disease_type"),
        "primary_site": record.get("primary_site"),
        "gender": safe_get(demographic, "gender"),
        "race": safe_get(demographic, "race"),
        "ethnicity": safe_get(demographic, "ethnicity"),
        "vital_status": vital_status,
        "days_to_birth": safe_get(demographic, "days_to_birth"),
        "age_at_index": safe_get(demographic, "age_at_index"),
        "days_to_death": days_to_death,
        "days_to_last_follow_up": days_to_last_follow_up,
        "survival_time": survival_time,
        "survival_event": survival_event,
        "tumor_stage": safe_get(diagnosis, "tumor_stage"),
        "tumor_grade": safe_get(diagnosis, "tumor_grade"),
        "primary_diagnosis": safe_get(diagnosis, "primary_diagnosis"),
        "ajcc_pathologic_stage": safe_get(diagnosis, "ajcc_pathologic_stage"),
        "ajcc_pathologic_t": safe_get(diagnosis, "ajcc_pathologic_t"),
        "ajcc_pathologic_n": safe_get(diagnosis, "ajcc_pathologic_n"),
        "ajcc_pathologic_m": safe_get(diagnosis, "ajcc_pathologic_m"),
    }


def build_clinical_table(records: list[dict]) -> pd.DataFrame:
    parsed = [parse_clinical_record(record) for record in records]
    return pd.DataFrame(parsed)


def build_master_table(primary_samples: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    master = primary_samples.merge(
        clinical,
        on="patient_id",
        how="left",
        validate="one_to_one",
    )

    return master


def summarize(master: pd.DataFrame) -> None:
    print("\n Master Clinical Table Summary ")
    print(f"Rows: {len(master)}")
    print(f"Unique patients: {master['patient_id'].nunique()}")

    print("\nVital status:")
    print(master["vital_status"].value_counts(dropna=False))

    print("\nSurvival event:")
    print(master["survival_event"].value_counts(dropna=False))

    print("\nMissing survival_time:")
    print(master["survival_time"].isna().sum())

    print("\nClinical columns preview:")
    clinical_cols = [
        "patient_id",
        "sample_barcode",
        "filename",
        "gender",
        "vital_status",
        "survival_time",
        "survival_event",
        "ajcc_pathologic_stage",
    ]
    available_cols = [col for col in clinical_cols if col in master.columns]
    print(master[available_cols].head().to_string(index=False))


def save(master: pd.DataFrame) -> None:
    out_dir = get_path_from_root("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)

    out_file = out_dir / "tcga_blca_master_clinical_table.csv"
    master.to_csv(out_file, index=False)

    print(f"\n[OK] Saved: {out_file}")


def main() -> int:
    try:
        primary_samples = load_primary_samples()
        clinical_records = load_clinical_json()
        clinical = build_clinical_table(clinical_records)

        master = build_master_table(primary_samples, clinical)

        summarize(master)
        save(master)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())