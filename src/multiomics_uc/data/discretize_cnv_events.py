from __future__ import annotations

import sys
import pandas as pd

from multiomics_uc.paths import get_path_from_root


def load_cnv_matrix() -> pd.DataFrame:
    path = get_path_from_root("data/processed/tcga_blca_cnv_gene_level_matrix.csv")
    return pd.read_csv(path, index_col=0)


def discretize_events(matrix: pd.DataFrame) -> dict[str, pd.DataFrame]:
    events = {
        "gain_cn_ge_3": (matrix >= 3).astype(int),
        "amplification_cn_ge_5": (matrix >= 5).astype(int),
        "deletion_cn_le_1": (matrix <= 1).astype(int),
        "deep_deletion_cn_eq_0": (matrix == 0).astype(int),
    }

    return events


def compute_event_burden(events: dict[str, pd.DataFrame]) -> pd.DataFrame:
    burden = pd.DataFrame(index=next(iter(events.values())).index)

    for event_name, event_matrix in events.items():
        burden[f"{event_name}_burden"] = event_matrix.sum(axis=1)

    burden = burden.reset_index().rename(columns={"index": "patient_id"})

    return burden


def summarize(events: dict[str, pd.DataFrame], burden: pd.DataFrame) -> None:
    print("\nCNV Event Discretization Summary")

    for event_name, event_matrix in events.items():
        total_events = int(event_matrix.sum().sum())
        mean_per_patient = event_matrix.sum(axis=1).mean()
        patients_with_event = (event_matrix.sum(axis=1) > 0).sum()

        print(f"\n{event_name}")
        print(f"Total gene-level events: {total_events}")
        print(f"Patients with at least one event: {patients_with_event}")
        print(f"Mean events per patient: {mean_per_patient:.2f}")

    print("\nBurden preview:")
    print(burden.head().to_string(index=False))


def save_outputs(events: dict[str, pd.DataFrame], burden: pd.DataFrame) -> None:
    processed_dir = get_path_from_root("data/processed")
    processed_dir.mkdir(parents=True, exist_ok=True)

    for event_name, event_matrix in events.items():
        out_file = processed_dir / f"tcga_blca_cnv_{event_name}_matrix.csv"
        event_matrix.to_csv(out_file)
        print(f"[OK] Saved {event_name} matrix: {out_file}")

    burden_file = processed_dir / "tcga_blca_cnv_event_burden.csv"
    burden.to_csv(burden_file, index=False)
    print(f"[OK] Saved CNV event burden: {burden_file}")


def main() -> int:
    try:
        matrix = load_cnv_matrix()

        events = discretize_events(matrix)
        burden = compute_event_burden(events)

        summarize(events, burden)
        save_outputs(events, burden)

        return 0

    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())