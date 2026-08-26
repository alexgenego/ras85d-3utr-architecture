#!/usr/bin/env python3
"""
Reference implementation of the final 63-branch enrichment classification.

Input table must contain:
node_or_branch, n_events_insertion, n_events_deletion, branch_length

The script intentionally does not reconstruct ARPIP events. It reproduces
branch fractions, expected counts, exact upper-tail binomial p-values,
BH-FDR and regime assignment after the canonical 63-branch table exists.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests


def classify(input_path: Path, output_path: Path):
    df = pd.read_csv(input_path)
    if len(df) != 63:
        raise ValueError(f"Expected 63 branches, found {len(df)}")
    if not (df["branch_length"] > 0).all():
        raise ValueError("All retained branches must have positive length")
    if "S.lebanonensis" in set(df["node_or_branch"].astype(str)):
        raise ValueError("Outgroup branch must be excluded")
    if "SYN_INGROUP_STEM" in set(df["node_or_branch"].astype(str)):
        raise ValueError("Synthetic ingroup stem must be excluded")

    df["branch_fraction"] = df["branch_length"] / df["branch_length"].sum()

    for event in ["insertion", "deletion"]:
        count_col = f"n_events_{event}"
        total = int(df[count_col].sum())
        df[f"expected_{event}s"] = total * df["branch_fraction"]
        df[f"OE_{event}s"] = np.where(
            df[f"expected_{event}s"] > 0,
            df[count_col] / df[f"expected_{event}s"],
            np.nan,
        )
        p = binom.sf(
            df[count_col].astype(int) - 1,
            total,
            df["branch_fraction"],
        )
        df[f"p_{event}s"] = p
        df[f"FDR_{event}s"] = multipletests(p, method="fdr_bh")[1]

    sig_i = df["FDR_insertions"] < 0.05
    sig_d = df["FDR_deletions"] < 0.05
    df["indel_regime"] = np.select(
        [sig_i & sig_d, sig_i, sig_d],
        ["bidirectional_turnover", "insertion_enriched", "deletion_enriched"],
        default="no_significant_enrichment",
    )

    if (int(df["n_events_insertion"].sum()),
        int(df["n_events_deletion"].sum())) != (124, 699):
        raise ValueError("Expected 124 insertions and 699 deletions")

    expected_regimes = {
        "bidirectional_turnover": 1,
        "insertion_enriched": 3,
        "deletion_enriched": 16,
        "no_significant_enrichment": 43,
    }
    observed = df["indel_regime"].value_counts().to_dict()
    if observed != expected_regimes:
        raise ValueError(f"Unexpected regime counts: {observed}")

    df.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    classify(args.input, args.output)
