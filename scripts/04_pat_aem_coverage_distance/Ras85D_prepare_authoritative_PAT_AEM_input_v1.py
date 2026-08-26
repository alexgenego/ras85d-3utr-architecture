#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the authoritative PAT/AEM input from 29A_3_01_final_object_registry.csv.

The script:
1. retains only authoritative FINAL records;
2. validates the final SM5/SM6 counts;
3. writes a MAT-compatible architecture table;
4. writes a minimal 1..1315 alignment3 master;
5. emits exact QC and inventory files.

It does not run PAT or AEM by itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


SPECIES = ["D.melanogaster", "D.yakuba", "D.virilis"]
EXPECTED = {
    "D.melanogaster": {
        "ECM_BLOCK": 17,
        "UGUA": 11,
        "PAS": 9,
        "CLEAVAGE_SITE": 12,
        "DOWNSTREAM_ELEMENT": 7,
        "RNA_STRUCTURE": 364,
    },
    "D.yakuba": {
        "ECM_BLOCK": 17,
        "UGUA": 9,
        "PAS": 9,
        "CLEAVAGE_SITE": 8,
        "DOWNSTREAM_ELEMENT": 10,
        "RNA_STRUCTURE": 334,
    },
    "D.virilis": {
        "ECM_BLOCK": 15,
        "UGUA": 12,
        "PAS": 13,
        "CLEAVAGE_SITE": 9,
        "DOWNSTREAM_ELEMENT": 18,
        "RNA_STRUCTURE": 437,
    },
}
TYPE_MAP = {
    "ECM": "ECM_BLOCK",
    "UGUA": "UGUA",
    "PAS": "PAS",
    "CLEAVAGE": "CLEAVAGE_SITE",
    "DSE": "DOWNSTREAM_ELEMENT",
    "STRUCTURE_SEGMENT": "RNA_STRUCTURE",
}
TRUE_VALUES = {"true", "t", "1", "yes", "y"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna("").astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def normalize_species(value: object) -> str:
    text = str(value).strip().replace("D. ", "D.")
    aliases = {
        "Dmelanogaster": "D.melanogaster",
        "Dyakuba": "D.yakuba",
        "Dvirilis": "D.virilis",
    }
    return aliases.get(text, text)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--alignment-length", type=int, default=1315)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(args.registry, low_memory=False)

    required = {
        "registry_object_id", "harmonized_object_id", "species", "object_type",
        "feature_group", "feature_class", "feature_subclass",
        "start", "end", "alignment_start", "alignment_end",
        "authoritative", "harmonization_status",
        "structure_model", "ecm_id", "annotation_value",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError("Final registry lacks columns: " + ", ".join(missing))

    source["species"] = source["species"].map(normalize_species)
    source["authoritative"] = as_bool(source["authoritative"])
    source["harmonization_status"] = (
        source["harmonization_status"].fillna("").astype(str).str.strip().str.upper()
    )

    final = source[
        source["authoritative"]
        & source["harmonization_status"].eq("FINAL")
        & source["species"].isin(SPECIES)
        & source["object_type"].isin(TYPE_MAP)
    ].copy()

    final["alignment_start"] = pd.to_numeric(
        final["alignment_start"], errors="raise"
    ).astype(int)
    final["alignment_end"] = pd.to_numeric(
        final["alignment_end"], errors="raise"
    ).astype(int)

    invalid = (
        (final["alignment_start"] < 1)
        | (final["alignment_end"] < final["alignment_start"])
        | (final["alignment_end"] > args.alignment_length)
    )
    if invalid.any():
        cols = [
            "registry_object_id", "species", "object_type",
            "alignment_start", "alignment_end",
        ]
        raise ValueError(
            "Invalid alignment3 coordinates:\n"
            + final.loc[invalid, cols].head(30).to_string(index=False)
        )

    if final["registry_object_id"].isna().any():
        final["registry_object_id"] = final["registry_object_id"].fillna(
            final["harmonized_object_id"]
        )
    if final["registry_object_id"].duplicated().any():
        duplicated = final.loc[
            final["registry_object_id"].duplicated(False),
            ["registry_object_id", "species", "object_type"],
        ]
        raise ValueError(
            "Duplicate registry_object_id:\n"
            + duplicated.head(30).to_string(index=False)
        )

    final["object_type_PAT"] = final["object_type"].map(TYPE_MAP)

    inventory = (
        final.groupby(["species", "object_type_PAT"])
        .size()
        .unstack(fill_value=0)
        .reindex(index=SPECIES, fill_value=0)
    )
    for object_type in [
        "ECM_BLOCK", "UGUA", "PAS", "CLEAVAGE_SITE",
        "DOWNSTREAM_ELEMENT", "RNA_STRUCTURE",
    ]:
        if object_type not in inventory.columns:
            inventory[object_type] = 0
    inventory = inventory[
        [
            "ECM_BLOCK", "UGUA", "PAS", "CLEAVAGE_SITE",
            "DOWNSTREAM_ELEMENT", "RNA_STRUCTURE",
        ]
    ]

    mismatches = []
    for species in SPECIES:
        for object_type, expected in EXPECTED[species].items():
            observed = int(inventory.loc[species, object_type])
            if observed != expected:
                mismatches.append(
                    {
                        "species": species,
                        "object_type": object_type,
                        "expected": expected,
                        "observed": observed,
                    }
                )
    if mismatches:
        raise ValueError(
            "Authoritative inventory mismatch:\n"
            + pd.DataFrame(mismatches).to_string(index=False)
        )

    mat = pd.DataFrame(
        {
            "object_id": final["registry_object_id"].astype(str),
            "harmonized_object_id": final["harmonized_object_id"].astype(str),
            "object_type": final["object_type_PAT"],
            "object_subtype": final["feature_subclass"].fillna(
                final["annotation_value"]
            ),
            "object_name": final["registry_object_id"].astype(str),
            "species": final["species"],
            "motif_family": final["ecm_id"].fillna(""),
            "sequence_start": pd.to_numeric(final["start"], errors="coerce"),
            "sequence_end": pd.to_numeric(final["end"], errors="coerce"),
            "alignment3_start": final["alignment_start"],
            "alignment3_end": final["alignment_end"],
            "coordinate_type": final["object_type_PAT"].map(
                lambda value: "boundary"
                if value == "CLEAVAGE_SITE"
                else "interval"
            ),
            "mapping_status": "authoritative_final",
            "structure_model": final["structure_model"].fillna(""),
            "structure_segment_id": final.get(
                "structure_segment_id", pd.Series("", index=final.index)
            ).fillna(""),
            "feature_group": final["feature_group"].fillna(""),
            "feature_class": final["feature_class"].fillna(""),
            "feature_subclass": final["feature_subclass"].fillna(""),
            "source": "29A_3_01_final_object_registry.csv",
            "confidence": "authoritative_final",
        }
    ).sort_values(
        ["species", "object_type", "alignment3_start", "alignment3_end", "object_id"]
    ).reset_index(drop=True)

    expected_mat_rows = 49 + 127 + 1135
    if len(mat) != expected_mat_rows:
        raise ValueError(
            f"Expected {expected_mat_rows} PAT/AEM objects, found {len(mat)}."
        )

    mat_path = args.output_dir / "Ras85D_Master_Architecture_Table_authoritative_v1.csv"
    compatibility_path = args.output_dir / "Ras85D_Master_Architecture_Table_v1_1.csv"
    master_path = args.output_dir / "Ras85D_alignment3_master_base.csv"
    inventory_path = args.output_dir / "Ras85D_authoritative_PAT_inventory_v1.csv"
    filtered_registry_path = (
        args.output_dir / "Ras85D_authoritative_PAT_AEM_registry_subset_v1.csv"
    )

    mat.to_csv(mat_path, index=False, encoding="utf-8-sig")
    mat.to_csv(compatibility_path, index=False, encoding="utf-8-sig")
    pd.DataFrame(
        {"alignment3_column": range(1, args.alignment_length + 1)}
    ).to_csv(master_path, index=False)
    inventory.reset_index().to_csv(inventory_path, index=False)
    final.to_csv(filtered_registry_path, index=False, encoding="utf-8-sig")

    qc = {
        "status": "PASS",
        "registry": str(args.registry),
        "registry_sha256": sha256(args.registry),
        "source_rows": int(len(source)),
        "authoritative_PAT_AEM_rows": int(len(final)),
        "MAT_rows": int(len(mat)),
        "alignment3_length": int(args.alignment_length),
        "expected_totals": {
            "ECM_BLOCK": 49,
            "APA": 127,
            "RNA_STRUCTURE": 1135,
            "all_PAT_AEM": 1311,
        },
        "observed_totals": {
            "ECM_BLOCK": int((mat["object_type"] == "ECM_BLOCK").sum()),
            "APA": int(
                mat["object_type"].isin(
                    ["UGUA", "PAS", "CLEAVAGE_SITE", "DOWNSTREAM_ELEMENT"]
                ).sum()
            ),
            "RNA_STRUCTURE": int((mat["object_type"] == "RNA_STRUCTURE").sum()),
            "all_PAT_AEM": int(len(mat)),
        },
        "inventory": inventory.reset_index().to_dict(orient="records"),
        "outputs": {
            "MAT": str(mat_path),
            "MAT_compatibility": str(compatibility_path),
            "alignment3_master": str(master_path),
            "inventory": str(inventory_path),
            "filtered_registry": str(filtered_registry_path),
        },
    }
    qc_path = args.output_dir / "Ras85D_authoritative_PAT_AEM_input_QC_v1.json"
    qc_path.write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 78)
    print("AUTHORITATIVE PAT/AEM INPUT — PASS")
    print("=" * 78)
    print(inventory.to_string())
    print()
    print("ECM:", qc["observed_totals"]["ECM_BLOCK"])
    print("APA:", qc["observed_totals"]["APA"])
    print("RNA structure:", qc["observed_totals"]["RNA_STRUCTURE"])
    print("All PAT/AEM objects:", qc["observed_totals"]["all_PAT_AEM"])
    print("QC:", qc_path)


if __name__ == "__main__":
    main()
