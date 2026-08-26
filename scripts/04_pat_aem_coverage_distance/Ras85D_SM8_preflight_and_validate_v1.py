#!/usr/bin/env python3
"""
Preflight the final Ras85D SM8 architecture input and validate transferred
Coverage Enrichment / Functional Distance outputs.

This wrapper does not replace the two final analysis modules. It prevents
silent use of legacy/non-final objects and checks the output schemas/QC values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED_ARCH = {
    "architecture_id", "feature_group", "feature_class",
    "feature_subclass", "architecture_label", "start", "end",
}
COVERAGE_REQUIRED = {
    "species", "domain", "test_family", "feature",
    "tested_region", "control_region",
    "perm_p_two_sided", "FDR_BH",
}
DISTANCE_REQUIRED = {
    "species", "feature", "landmark",
    "observed_median_distance_nt", "null_median_distance_nt",
    "perm_p_two_sided", "FDR_BH",
}


def read_table(path: Path, sheet: str | None = None) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet or 0)
    raise ValueError(f"Unsupported table: {path}")


def bool_series(series: pd.Series) -> pd.Series:
    return series.fillna(False).astype(str).str.strip().str.lower().isin(
        {"true", "1", "yes", "y", "pass", "authoritative", "final"}
    )


def preflight_architecture(path: Path, utr_end: int) -> pd.DataFrame:
    df = read_table(path)
    missing = sorted(REQUIRED_ARCH - set(df.columns))
    if missing:
        raise ValueError("Architecture missing columns: " + ", ".join(missing))

    if "species" not in df.columns:
        df["species"] = "D.melanogaster"

    df["architecture_id"] = df["architecture_id"].astype(str).str.strip()
    if df["architecture_id"].eq("").any() or df["architecture_id"].duplicated().any():
        raise ValueError("architecture_id must be non-empty and unique")

    df["start"] = pd.to_numeric(df["start"], errors="raise").astype(int)
    df["end"] = pd.to_numeric(df["end"], errors="raise").astype(int)

    invalid = (
        (df["start"] < 1)
        | (df["end"] < df["start"])
        | (df["end"] > utr_end)
    )
    if invalid.any():
        raise ValueError(
            "Intervals outside 1-based inclusive bounds:\n"
            + df.loc[invalid, ["architecture_id", "start", "end"]]
            .head(20).to_string(index=False)
        )

    # Detect common authority/status columns. If present, no non-final row may remain.
    authority_columns = [
        c for c in df.columns
        if c.lower() in {
            "authoritative", "is_authoritative", "final_object",
            "is_final", "analysis_ready", "include_in_analysis",
        }
    ]
    status_columns = [
        c for c in df.columns
        if c.lower() in {
            "status", "record_status", "object_status", "harmonization_status",
        }
    ]

    for column in authority_columns:
        bad = ~bool_series(df[column])
        if bad.any():
            raise ValueError(
                f"Non-authoritative rows remain according to {column}: "
                f"{int(bad.sum())}"
            )

    for column in status_columns:
        text = df[column].fillna("").astype(str).str.upper()
        bad = text.str.contains("LEGACY|NON_FINAL|NON-FINAL|EXCLUDED", regex=True)
        if bad.any():
            raise ValueError(
                f"Legacy/non-final rows remain according to {column}: "
                f"{int(bad.sum())}"
            )

    print("Architecture preflight PASS")
    print("Rows:", len(df))
    print("Species:", sorted(df["species"].astype(str).unique()))
    print("Coordinate range:", int(df["start"].min()), int(df["end"].max()))
    print("Explicit UTR end:", utr_end)
    return df


def validate_qc(path: Path, expected_name: str) -> dict:
    qc = json.loads(path.read_text(encoding="utf-8"))
    if expected_name.lower() not in str(qc.get("test", "")).lower():
        raise ValueError(f"Unexpected QC test label in {path}")
    if int(qc.get("n_permutations", -1)) != 10000:
        raise ValueError(f"Expected 10000 permutations in {path}")
    if int(qc.get("seed", -1)) != 20260718:
        raise ValueError(f"Expected seed 20260718 in {path}")
    return qc


def validate_coverage(results_path: Path, qc_path: Path) -> None:
    df = read_table(results_path, sheet="All_results")
    missing = sorted(COVERAGE_REQUIRED - set(df.columns))
    if missing:
        raise ValueError("Coverage results missing: " + ", ".join(missing))
    if "analysis_type" in df.columns and df["analysis_type"].astype(str).str.contains(
        "DISTANCE", case=False, regex=False
    ).any():
        raise ValueError("Coverage-only output contains distance tests")
    if not ((pd.to_numeric(df["FDR_BH"], errors="coerce") >= 0)
            & (pd.to_numeric(df["FDR_BH"], errors="coerce") <= 1)).all():
        raise ValueError("Coverage FDR values outside [0,1]")

    qc = validate_qc(qc_path, "Coverage Enrichment Test")
    if int(qc.get("tests", -1)) != len(df):
        raise ValueError("Coverage QC test count does not match result rows")
    if int(qc.get("apa_window_nt", -1)) != 25:
        raise ValueError("Coverage APA window is not 25 nt")

    print("Coverage output PASS")
    print("Rows:", len(df))
    print("Nominal p<0.05:", int((pd.to_numeric(df["perm_p_two_sided"]) < 0.05).sum()))
    print("FDR<0.05:", int((pd.to_numeric(df["FDR_BH"]) < 0.05).sum()))


def validate_distance(results_path: Path, qc_path: Path) -> None:
    df = read_table(results_path, sheet="All_results")
    missing = sorted(DISTANCE_REQUIRED - set(df.columns))
    if missing:
        raise ValueError("Distance results missing: " + ", ".join(missing))

    landmarks = set(df["landmark"].astype(str))
    if not landmarks.issubset({"PAS", "CLEAVAGE"}):
        raise ValueError(f"Unexpected primary landmarks: {sorted(landmarks)}")
    self_rows = (
        (df["feature"].astype(str) == "PAS") & (df["landmark"].astype(str) == "PAS")
    ) | (
        (df["feature"].astype(str) == "CLEAVAGE")
        & (df["landmark"].astype(str) == "CLEAVAGE")
    )
    if self_rows.any():
        raise ValueError("Primary distance output contains self-comparisons")

    qc = validate_qc(qc_path, "Functional Distance Test")
    if int(qc.get("tests", -1)) != len(df):
        raise ValueError("Distance QC test count does not match result rows")
    if not bool(qc.get("exclude_self", False)):
        raise ValueError("Distance QC does not confirm exclude_self=True")

    print("Distance output PASS")
    print("Rows:", len(df))
    print("Nominal p<0.05:", int((pd.to_numeric(df["perm_p_two_sided"]) < 0.05).sum()))
    print("FDR<0.05:", int((pd.to_numeric(df["FDR_BH"]) < 0.05).sum()))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--architecture", type=Path)
    parser.add_argument("--utr-end", type=int, default=964)
    parser.add_argument("--coverage-results", type=Path)
    parser.add_argument("--coverage-qc", type=Path)
    parser.add_argument("--distance-results", type=Path)
    parser.add_argument("--distance-qc", type=Path)
    args = parser.parse_args()

    if args.architecture:
        preflight_architecture(args.architecture, args.utr_end)
    if args.coverage_results or args.coverage_qc:
        if not (args.coverage_results and args.coverage_qc):
            raise ValueError("Provide both --coverage-results and --coverage-qc")
        validate_coverage(args.coverage_results, args.coverage_qc)
    if args.distance_results or args.distance_qc:
        if not (args.distance_results and args.distance_qc):
            raise ValueError("Provide both --distance-results and --distance-qc")
        validate_distance(args.distance_results, args.distance_qc)

    if not any([
        args.architecture,
        args.coverage_results,
        args.distance_results,
    ]):
        parser.error("No validation target supplied")


if __name__ == "__main__":
    main()
