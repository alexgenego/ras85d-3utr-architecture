#!/usr/bin/env python3
"""Ras85D multiloop-state dispersion analysis.

Inputs:
  1. 29A_3_08B_Tissue_Resolved_APA_Cleavage_Efficiency_v2.xlsx
  2. 29A_3_09B_multiloop_state_data.xlsx

The script calculates site-centered empirical-logit dispersion for five
tissue/developmental contexts and tests structural-state weighted variances
using within-species permutations of complete five-state profiles.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd


STATES = [
    "SAME_MULTILOOP_INTERVAL",
    "DIFFERENT_INTERVALS_SAME_MULTILOOP",
    "MULTILOOP_TO_ARM",
    "SAME_MULTILOOP_ARM",
    "DIFFERENT_MULTILOOP_ARMS",
]
CONTEXTS = [
    "HEAD",
    "OVARY",
    "TESTIS",
    "EARLY_EMBRYO_0_12H",
    "LATE_EMBRYO_12_24H",
]
PREFIX = "consensus_NUCLEOTIDE_FRACTIONAL_ALL_PAIRS_fraction_"


def empirical_logit(k: float, n: float) -> float:
    return math.log((k + 0.5) / (n - k + 0.5))


def weighted_variance(centered: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights[:, None] * centered**2) /
                 (centered.shape[1] * np.sum(weights)))


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    order = np.argsort(pvalues)
    ranked = pvalues[order]
    q = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--usage", required=True, type=Path)
    parser.add_argument("--structure", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    usage = pd.read_excel(
        args.usage, sheet_name="Counts_and_Fractions", header=3
    )
    structure = pd.read_excel(
        args.structure, sheet_name="Site_Metrics", header=3
    )

    columns = ["cleavage_id", "species", "cleavage_position"] + [
        PREFIX + state for state in STATES
    ]
    structure = structure[columns].copy()
    merged = usage.merge(structure, on=["cleavage_id", "species"], how="inner")

    site_rows = []
    for cleavage_id, frame in merged.groupby("cleavage_id", sort=False):
        by_context = frame.set_index("tissue_group")
        if not set(CONTEXTS).issubset(by_context.index):
            continue
        state_values = by_context.iloc[0][[PREFIX + s for s in STATES]].astype(float)
        if not np.isfinite(state_values).all():
            continue

        logits = np.array([
            empirical_logit(
                float(by_context.loc[c, "count"]),
                float(by_context.loc[c, "mapped_sample_total"]),
            )
            for c in CONTEXTS
        ])
        site_rows.append({
            "cleavage_id": cleavage_id,
            "species": by_context.iloc[0]["species"],
            "cleavage_position": by_context.iloc[0]["cleavage_position"],
            "logits": logits,
            "states": state_values.to_numpy(float),
        })

    species = np.array([r["species"] for r in site_rows], dtype=object)
    logits = np.vstack([r["logits"] for r in site_rows])
    state_matrix = np.vstack([r["states"] for r in site_rows])
    centered = logits - logits.mean(axis=1, keepdims=True)

    observed = np.array([
        weighted_variance(centered, state_matrix[:, j])
        for j in range(len(STATES))
    ])
    rng = np.random.default_rng(args.seed)
    species_indices = [
        np.where(species == sp)[0] for sp in np.unique(species)
    ]
    null = np.empty((args.permutations, len(STATES)))
    contrast = np.empty(args.permutations)

    for b in range(args.permutations):
        permuted = state_matrix.copy()
        for idx in species_indices:
            permuted[idx] = state_matrix[rng.permutation(idx)]
        null[b] = [
            weighted_variance(centered, permuted[:, j])
            for j in range(len(STATES))
        ]
        contrast[b] = null[b, 4] - null[b, 0]

    p_high = (np.sum(null >= observed, axis=0) + 1) / (args.permutations + 1)
    p_low = (np.sum(null <= observed, axis=0) + 1) / (args.permutations + 1)
    p_two = np.minimum(1.0, 2 * np.minimum(p_high, p_low))

    results = pd.DataFrame({
        "state": STATES,
        "observed_weighted_variance": observed,
        "null_mean": null.mean(axis=0),
        "observed_null_ratio": observed / null.mean(axis=0),
        "p_high": p_high,
        "p_low": p_low,
        "p_two": p_two,
        "BH_FDR_two_sided": bh_fdr(p_two),
    })
    results.to_csv(args.outdir / "primary_five_states.csv", index=False)

    observed_contrast = observed[4] - observed[0]
    contrast_p = (
        np.sum(contrast >= observed_contrast) + 1
    ) / (args.permutations + 1)
    pd.DataFrame([{
        "contrast": "DIFFERENT_MULTILOOP_ARMS_minus_SAME_MULTILOOP_INTERVAL",
        "observed_difference": observed_contrast,
        "permutation_p_high": contrast_p,
    }]).to_csv(args.outdir / "primary_contrast.csv", index=False)

    np.savez_compressed(
        args.outdir / "null_distributions.npz",
        states=np.array(STATES),
        null_state_variances=null,
        null_primary_contrast=contrast,
    )


if __name__ == "__main__":
    main()
