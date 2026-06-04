from __future__ import annotations

import sys
import pandas as pd
from scipy.stats import percentileofscore

from multiomics_uc.paths import get_path_from_root


PROGRAMS = {
    "PC1": "immune_infiltration",
    "PC2": "angiogenic_stroma",
    "PC3": "proliferation",
    "PC4": "motility_cytoskeleton",
    "PC6": "immune_activation",
    "PC8": "invasion_signaling",
    "PC9": "differentiation",
}


def load_data() -> pd.DataFrame:
    path = get_path_from_root(
        "data/processed/integrated_patient_feature_table_clean.csv"
    )
    return pd.read_csv(path)


def compute_percentiles(df: pd.DataFrame) -> pd.DataFrame:

    result = pd.DataFrame()
    result["patient_id"] = df["patient_id"]

    for pc, label in PROGRAMS.items():

        values = df[pc]

        result[f"{label}_percentile"] = [
            percentileofscore(values, x)
            for x in values
        ]

    # -------------------------
    # Mutation burden percentile
    # -------------------------

    if "mutation_burden" in df.columns:

        values = df["mutation_burden"].dropna()

        result["mutation_burden_percentile"] = [
            percentileofscore(values, x)
            if pd.notna(x)
            else None
            for x in df["mutation_burden"]
        ]

    return result

def main() -> int:

    try:

        df = load_data()

        percentiles = compute_percentiles(df)

        out_path = get_path_from_root(
            "data/processed/program_percentiles.csv"
        )

        percentiles.to_csv(
            out_path,
            index=False,
        )

        print("\nProgram Percentiles")
        print(percentiles.head())

        print(f"\n[OK] Saved: {out_path}")

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())