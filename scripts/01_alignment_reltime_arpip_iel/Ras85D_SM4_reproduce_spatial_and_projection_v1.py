#!/usr/bin/env python3
"""
Reproduce the Ras85D SM4 spatial statistics from the bundled 34-IEL CSV.

Inputs:
  Ras85D_SM4_IEL_regions_34_minimal_v1.csv
  Ras85D_SM4_Dmel_loci_27_v1.csv

The Monte Carlo procedure samples 34 unique integer centers from 1..2247.
Two edge spacings are included: 0.5 to the first center and the last center
to 2247.5. The CV is the sample standard deviation of the 35 spacings
divided by their mean. Upper-tail p-values use the +1 correction.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


def read_iel(path: Path):
    with path.open(encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    centers = np.array([float(row["alignment37_center"]) for row in rows])
    lengths = np.array([int(row["alignment37_length"]) for row in rows])
    return rows, centers, lengths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--iel", type=Path, required=True)
    parser.add_argument("--n-permutations", type=int, default=100000)
    parser.add_argument("--seed", type=int, default=20260716)
    args = parser.parse_args()

    rows, centers, lengths = read_iel(args.iel)
    if len(rows) != 34:
        raise ValueError(f"Expected 34 IEL, found {len(rows)}")
    if int(lengths.sum()) != 123:
        raise ValueError(f"Expected total IEL span 123, found {lengths.sum()}")

    centers = np.sort(centers)
    spacings = np.diff(np.concatenate(([0.5], centers, [2247.5])))
    observed_cv = spacings.std(ddof=1) / spacings.mean()
    observed_max = spacings.max()

    rng = np.random.default_rng(args.seed)
    ge_cv = 0
    ge_max = 0

    for _ in range(args.n_permutations):
        random_centers = np.sort(
            rng.choice(np.arange(1, 2248), size=34, replace=False).astype(float)
        )
        random_spacings = np.diff(
            np.concatenate(([0.5], random_centers, [2247.5]))
        )
        cv = random_spacings.std(ddof=1) / random_spacings.mean()
        maximum = random_spacings.max()
        ge_cv += int(cv >= observed_cv)
        ge_max += int(maximum >= observed_max)

    p_cv = (ge_cv + 1) / (args.n_permutations + 1)
    p_max = (ge_max + 1) / (args.n_permutations + 1)

    print("n_IEL:", len(rows))
    print("total_span:", int(lengths.sum()))
    print("median_length:", float(np.median(lengths)))
    print("length_range:", int(lengths.min()), int(lengths.max()))
    print("spacing_CV:", observed_cv)
    print("spacing_CV_p:", p_cv)
    print("maximum_spacing:", observed_max)
    print("maximum_spacing_p:", p_max)


if __name__ == "__main__":
    main()
