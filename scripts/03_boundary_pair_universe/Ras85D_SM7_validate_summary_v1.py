#!/usr/bin/env python3
"""Validate the compact Ras85D SM7 summary tables."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def numeric_map(rows, key, value):
    return {row[key]: float(row[value]) for row in rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bpu", type=Path, required=True)
    parser.add_argument("--interfaces", type=Path, required=True)
    parser.add_argument("--states", type=Path, required=True)
    parser.add_argument("--type-change", type=Path, required=True)
    args = parser.parse_args()

    bpu = numeric_map(read_rows(args.bpu), "metric", "value")
    interfaces = numeric_map(read_rows(args.interfaces), "metric", "value")
    states = read_rows(args.states)
    type_change = read_rows(args.type_change)

    assert int(bpu["Pairwise-enabled QC-pass objects"]) == 525
    assert int(bpu["All unordered object pairs"]) == 137550
    assert int(bpu["Eligible pairs"]) == 107321
    assert int(bpu["Ineligible pairs"]) == 30229
    assert int(bpu["Observed interfaces mapped to BPU"]) == 871
    assert int(bpu["All unordered object pairs"]) == 525 * 524 // 2
    assert int(bpu["Eligible pairs"] + bpu["Ineligible pairs"]) == 137550

    assert int(interfaces["Observed strict functional interfaces"]) == 109
    assert int(interfaces["Homologous interface families"]) == 78
    assert int(interfaces["Species-state matrix rows"]) == 234
    assert int(interfaces["Exact type/configuration conserved"]) == 29
    assert int(interfaces["Singleton, not comparable"]) == 48
    assert int(interfaces["Type/configuration change"]) == 1
    assert int(interfaces["Critical QC failures"]) == 0

    assert len(states) == 5
    assert len(type_change) == 1
    assert int(float(type_change[0]["minimum_type_changes"])) == 1
    assert int(float(type_change[0]["equally_parsimonious_placements"])) == 3

    print("SM7 summary QC passed")
    print("BPU:", 525, "objects;", 137550, "pairs;", 107321, "eligible")
    print("Observed BPU interfaces:", 871)
    print("Functional interfaces:", 109)
    print("Homologous families:", 78)
    print("Type changes:", 1)


if __name__ == "__main__":
    main()
