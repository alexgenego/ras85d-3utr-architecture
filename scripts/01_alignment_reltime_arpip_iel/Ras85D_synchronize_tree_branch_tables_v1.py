#!/usr/bin/env python3
"""
Synchronize Ras85D tree names, legacy V-node IDs and table branch fields.

The script does not recompute branch statistics. For final statistics it treats
Table_S6_Branch_level_statistics.xlsx ("63 branches") as the authoritative
source, because expected counts, O/E and FDR must be recalculated after changing
the branch universe.

Usage in Colab:
    python Ras85D_synchronize_tree_branch_tables_v1.py \
        --crosswalk Ras85D_branch_ID_crosswalk_v1.csv \
        --input Table_S6_Branch_level_statistics.xlsx \
        --output Table_S6_Branch_level_statistics_canonical.xlsx

It also works with Table_1_Significant_indel_regimes.xlsx,
Table_S5_IEL_regions_full.xlsx and CSV/TSV files.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

NAME_CORRECTIONS = {
    "D.willistony": "D.willistoni",
    "D.montana85": "D.montana",
    "D.kanecoi": "D.kanekoi",
    "D.littoralis74": "D.littoralis",
    "D.lummmei": "D.lummei",
}

BRANCH_COLUMN_NAMES = {
    "node_or_branch", "Branch", "branch", "branch_id",
    "deletion_branches", "insertion_branches",
}

DESCENDANT_COLUMN_NAMES = {
    "descendant_taxa", "Species", "species", "taxon", "Taxon",
}

def load_crosswalk(path: Path) -> dict[str, str]:
    cw = pd.read_csv(path)
    return dict(zip(cw["legacy_branch_id"], cw["canonical_branch_id"]))

def replace_tokens(value, mapping):
    if pd.isna(value):
        return value
    text = str(value)

    for old, new in NAME_CORRECTIONS.items():
        text = re.sub(rf"(?<![A-Za-z0-9_.]){re.escape(old)}(?![A-Za-z0-9_.])", new, text)

    # Longest tokens first prevents partial replacements.
    for old in sorted(mapping, key=len, reverse=True):
        new = mapping[old]
        text = re.sub(rf"(?<![A-Za-z0-9_.]){re.escape(old)}(?![A-Za-z0-9_.])", new, text)

    return text

def synchronize_frame(df: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    out = df.copy()

    for col in out.columns:
        if (
            col in BRANCH_COLUMN_NAMES
            or col in DESCENDANT_COLUMN_NAMES
            or "branch" in str(col).lower()
            or "tax" in str(col).lower()
            or "species" in str(col).lower()
        ):
            out[col] = out[col].map(lambda x: replace_tokens(x, mapping))

    # Add review flags without discarding legacy evidence.
    text_cols = [c for c in out.columns if out[c].dtype == object]
    if text_cols:
        joined = out[text_cols].fillna("").astype(str).agg(" ; ".join, axis=1)
        out["contains_synthetic_ingroup_stem"] = joined.str.contains(
            "SYN_INGROUP_STEM", regex=False
        )
        out["contains_root_reference"] = joined.str.contains(
            r"(?<![A-Za-z0-9_])ROOT(?![A-Za-z0-9_])", regex=True
        )

    return out

def validate_final_63(df: pd.DataFrame) -> None:
    if len(df) != 63:
        raise ValueError(f"Expected 63 branches, found {len(df)}.")

    if "branch_length" in df:
        if not (pd.to_numeric(df["branch_length"], errors="coerce") > 0).all():
            raise ValueError("The final 63-branch table contains a non-positive branch length.")

    if "node_or_branch" in df:
        forbidden = {"S.lebanonensis", "SYN_INGROUP_STEM"}
        found = forbidden.intersection(set(df["node_or_branch"].astype(str)))
        if found:
            raise ValueError(f"Forbidden branches in final 63 universe: {sorted(found)}")

    if {"n_events_insertion", "n_events_deletion"}.issubset(df.columns):
        ins = int(pd.to_numeric(df["n_events_insertion"]).sum())
        dele = int(pd.to_numeric(df["n_events_deletion"]).sum())
        if (ins, dele) != (124, 699):
            raise ValueError(
                f"Expected event totals 124 insertions / 699 deletions; found {ins} / {dele}."
            )

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--crosswalk", type=Path, required=True)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    mapping = load_crosswalk(args.crosswalk)
    suffix = args.input.suffix.lower()

    if suffix in {".csv", ".tsv"}:
        sep = "\t" if suffix == ".tsv" else ","
        df = pd.read_csv(args.input, sep=sep)
        out = synchronize_frame(df, mapping)
        out.to_csv(args.output, sep=sep, index=False)
        return

    if suffix not in {".xlsx", ".xlsm"}:
        raise ValueError(f"Unsupported input: {suffix}")

    sheets = pd.read_excel(args.input, sheet_name=None)
    synced = {}

    for sheet_name, df in sheets.items():
        out = synchronize_frame(df, mapping)

        if sheet_name.strip().lower() == "63 branches":
            validate_final_63(out)

        synced[sheet_name] = out

    with pd.ExcelWriter(args.output, engine="openpyxl") as writer:
        for sheet_name, df in synced.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)

    print(f"Saved: {args.output}")

if __name__ == "__main__":
    main()
