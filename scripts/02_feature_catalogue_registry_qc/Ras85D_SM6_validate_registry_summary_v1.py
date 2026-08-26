#!/usr/bin/env python3
"""
Validate the reviewer-facing SM6 summary tables.

This validates count consistency. The complete harmonized 1554-row source
registry remains in the authoritative 29A.3 workbooks and is not reconstructed
from summary values.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", type=Path, required=True)
    parser.add_argument("--types", type=Path, required=True)
    parser.add_argument("--scope", type=Path, required=True)
    args = parser.parse_args()

    species = read_csv(args.species)
    types = read_csv(args.types)
    scope = read_csv(args.scope)

    functional = sum(int(row["primary_functional_objects"]) for row in species)
    structures = sum(int(row["structure_segments"]) for row in species)
    authoritative = sum(int(row["authoritative_objects_total"]) for row in species)

    type_objects = sum(int(row["authoritative_objects"]) for row in types)
    refined = sum(
        int(row["refined_orthology_groups"])
        for row in types
        if row["refined_orthology_groups"].strip()
    )
    scope_groups = sum(int(row["refined_groups"]) for row in scope)

    assert functional == 390
    assert structures == 1135
    assert authoritative == 1525
    assert type_objects == 1525
    assert refined == 214
    assert scope_groups == 214
    assert len(species) == 3

    print("SM6 summary QC passed")
    print("Primary functional objects:", functional)
    print("Structure segments:", structures)
    print("Authoritative objects:", authoritative)
    print("Refined orthology groups:", refined)
    print("Species:", ", ".join(row["species"] for row in species))


if __name__ == "__main__":
    main()
