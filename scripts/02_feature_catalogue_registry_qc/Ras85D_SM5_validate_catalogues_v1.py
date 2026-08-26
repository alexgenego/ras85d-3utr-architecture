#!/usr/bin/env python3
"""
Validate the reviewer-facing Ras85D SM5 catalogues.

This script validates the assembled Table S3 evidence. It does not reproduce
the original CENSOR/Repbase or MEME/MAST searches, because the raw discovery
directories are not included in the current available package.
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path


SPECIES_LENGTHS = {
    "D. melanogaster": 964,
    "D. yakuba": 968,
    "D. virilis": 1150,
}


def read_csv(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mobile-elements", type=Path, required=True)
    parser.add_argument("--microsatellites", type=Path, required=True)
    parser.add_argument("--ecm-groups", type=Path, required=True)
    parser.add_argument("--ecm-occurrences", type=Path, required=True)
    args = parser.parse_args()

    te = read_csv(args.mobile_elements)
    ms = read_csv(args.microsatellites)
    groups = read_csv(args.ecm_groups)
    occurrences = read_csv(args.ecm_occurrences)

    assert len(te) == 3
    assert {row["status"] for row in te} == {"accepted candidate trace"}
    assert len(groups) == 17
    assert len(occurrences) == 49

    counts = {species: 0 for species in SPECIES_LENGTHS}
    reverse = []

    for row in occurrences:
        species = row["species"]
        start = int(row["start"])
        end = int(row["end"])
        length = int(row["length"])

        assert 1 <= start <= end <= SPECIES_LENGTHS[species]
        assert length == end - start + 1
        counts[species] += 1

        if row["orientation"] == "reverse_complement":
            reverse.append((row["ECM_group"], species))

    assert counts == {
        "D. melanogaster": 17,
        "D. yakuba": 17,
        "D. virilis": 15,
    }
    assert reverse == [("ECM_11", "D. virilis")]

    dvir_groups = {
        row["ECM_group"]
        for row in occurrences
        if row["species"] == "D. virilis"
    }
    assert "ECM_16" not in dvir_groups
    assert "ECM_17" not in dvir_groups

    assert any(
        row["journal_status"] == "strong comparative example"
        for row in ms
    )
    assert all(
        "pending" in row["coordinate_status"]
        or "documented" in row["coordinate_status"]
        for row in ms
    )

    print("SM5 catalogue QC passed")
    print("TE candidates:", len(te))
    print("Microsatellite evidence rows:", len(ms))
    print("ECM groups:", len(groups))
    print("ECM occurrences:", len(occurrences))
    print("Species counts:", counts)
    print("Reverse-complement occurrence:", reverse[0])


if __name__ == "__main__":
    main()
