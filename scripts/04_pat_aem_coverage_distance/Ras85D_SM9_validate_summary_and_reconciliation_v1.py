#!/usr/bin/env python3
"""Validate compact SM9 evidence and flag unresolved PAT/AEM input mismatches."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pat-edges", type=Path, required=True)
    p.add_argument("--pat-inventory", type=Path, required=True)
    p.add_argument("--iel-ecm", type=Path, required=True)
    p.add_argument("--reconciliation", type=Path, required=True)
    args = p.parse_args()

    edges = read(args.pat_edges)
    inventory = read(args.pat_inventory)
    tests = read(args.iel_ecm)
    rec = read(args.reconciliation)

    assert len(edges) == 3
    significant = [row for row in edges if row["statistical_class"] == "significant_enrichment"]
    assert len(significant) == 1
    assert significant[0]["edge_id"] == "RNA_STRUCTURE__APA"
    assert abs(float(significant[0]["effect_ratio"]) - 5.282327) < 1e-6
    assert abs(float(significant[0]["q_enrichment"]) - 0.001521) < 1e-9

    assert len(inventory) == 3
    assert sum(int(row["PAT_ECM"]) for row in inventory) == 49
    assert sum(int(row["PAT_structure_total"]) for row in inventory) == 1092

    assert len(tests) == 5
    assert all(float(row["FDR"]) < 0.05 for row in tests)

    assert len(rec) == 3
    structure_delta = sum(int(row["delta_structure"]) for row in rec)
    assert structure_delta == -43
    assert any(int(row["delta_APA"]) != 0 for row in rec)

    print("SM9 compact evidence QC PASS")
    print("PAT significant edges: 1 (RNA_STRUCTURE__APA)")
    print("PAT ECM occurrences:", 49)
    print("PAT structure segments:", 1092)
    print("Authoritative structure deficit:", 43)
    print("Blocking reconciliation remains: YES")


if __name__ == "__main__":
    main()
