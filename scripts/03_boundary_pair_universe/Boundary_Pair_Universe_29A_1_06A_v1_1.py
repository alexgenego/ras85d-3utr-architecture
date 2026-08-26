#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
29A.1.06A v1.1 — Boundary Pair Universe (BPU)

Builds the frozen universe of unordered pairs from the canonical Ras85D architecture.
The v1.1 logic uses the architecture table's own control columns:
  * pairwise_enabled — primary inclusion flag;
  * pairwise_category — canonical object category;
  * coordinate_valid, length_matches_coordinates, within_utr_bounds — QC flags.

Coordinates are 1-based inclusive.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

STAGE = "29A.1.06A"
VERSION = "1.1"
REQUIRED = {
    "architecture_id", "start", "end", "species",
    "pairwise_enabled", "pairwise_category",
    "coordinate_valid", "length_matches_coordinates", "within_utr_bounds",
}
TRUE_VALUES = {"true", "t", "1", "yes", "y"}


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)
    raise ValueError(f"Unsupported input format: {suffix}")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_observed_pairs(path: Optional[Path]) -> set[tuple[str, str]]:
    if path is None:
        return set()
    df = read_table(path)
    alternatives = [
        ("object_a", "object_b"),
        ("architecture_id_a", "architecture_id_b"),
        ("source", "target"),
        ("node_a", "node_b"),
    ]
    pair_cols = next((x for x in alternatives if set(x).issubset(df.columns)), None)
    if pair_cols is None:
        raise ValueError("Observed-pairs table lacks a supported pair-column combination.")
    pairs = set()
    for a, b in df[list(pair_cols)].itertuples(index=False, name=None):
        a, b = str(a).strip(), str(b).strip()
        if a and b and a != b:
            pairs.add(tuple(sorted((a, b))))
    return pairs


def geometry(a_start: int, a_end: int, b_start: int, b_end: int) -> dict:
    overlap = max(0, min(a_end, b_end) - max(a_start, b_start) + 1)
    if (a_start, a_end) <= (b_start, b_end):
        left_s, left_e, right_s, right_e = a_start, a_end, b_start, b_end
        order = "A_LEFT_OF_B"
    else:
        left_s, left_e, right_s, right_e = b_start, b_end, a_start, a_end
        order = "B_LEFT_OF_A"

    if overlap > 0:
        gap = 0
        edge_distance = 0
        relation = "OVERLAP"
    else:
        gap = max(0, right_s - left_e - 1)
        edge_distance = max(0, right_s - left_e)
        relation = "TOUCHING" if gap == 0 else "SEPARATED"

    exact = a_start == b_start and a_end == b_end
    if exact:
        containment = "COEXTENSIVE"
    elif a_start <= b_start and a_end >= b_end:
        containment = "A_CONTAINS_B"
    elif b_start <= a_start and b_end >= a_end:
        containment = "B_CONTAINS_A"
    elif overlap:
        containment = "PARTIAL_OVERLAP"
    else:
        containment = "DISJOINT"

    return {
        "geometric_order": order,
        "relation": relation,
        "containment": containment,
        "overlap_nt": int(overlap),
        "gap_nt": int(gap),
        "boundary_distance_nt": int(edge_distance),
        "span_start": int(min(a_start, b_start)),
        "span_end": int(max(a_end, b_end)),
        "span_length_nt": int(max(a_end, b_end) - min(a_start, b_start) + 1),
        "exact_coordinate_match": bool(exact),
    }


def prepare(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    missing = sorted(REQUIRED - set(df.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))

    x = df.copy()
    x["architecture_id"] = x["architecture_id"].astype(str).str.strip()
    if x["architecture_id"].eq("").any() or x["architecture_id"].duplicated().any():
        raise ValueError("architecture_id must be non-empty and unique.")

    x["start"] = pd.to_numeric(x["start"], errors="raise").astype(int)
    x["end"] = pd.to_numeric(x["end"], errors="raise").astype(int)
    x["pairwise_enabled"] = as_bool(x["pairwise_enabled"])
    for c in ["coordinate_valid", "length_matches_coordinates", "within_utr_bounds"]:
        x[c] = as_bool(x[c])
    x["pairwise_category"] = x["pairwise_category"].astype(str).str.strip()
    x["species"] = x["species"].astype(str).str.strip()
    x["structure_model"] = x.get("structure_model", "").fillna("").astype(str).str.strip().str.upper()
    x["feature_group"] = x.get("feature_group", "").fillna("").astype(str).str.strip()

    qc_ok = x[["coordinate_valid", "length_matches_coordinates", "within_utr_bounds"]].all(axis=1)
    x["bpu_input_status"] = np.select(
        [~x["pairwise_enabled"], x["pairwise_enabled"] & ~qc_ok],
        ["EXCLUDED_PAIRWISE_DISABLED", "EXCLUDED_QC_FAILURE"],
        default="INCLUDED",
    )
    included = x[x["bpu_input_status"].eq("INCLUDED")].copy()
    excluded = x[~x["bpu_input_status"].eq("INCLUDED")].copy()

    if included["pairwise_category"].eq("").any():
        raise ValueError("Included objects contain empty pairwise_category values.")
    if ((included["start"] < 1) | (included["end"] < included["start"])).any():
        raise ValueError("Included objects contain invalid numeric coordinates despite QC flags.")

    included["length_nt"] = included["end"] - included["start"] + 1
    included["center"] = (included["start"] + included["end"]) / 2.0
    included["architecture_layer"] = np.where(
        included["feature_group"].str.upper().eq("RNA_STRUCTURE"),
        "RNA_STRUCTURE",
        included["feature_group"].str.upper().replace("", "UNASSIGNED"),
    )
    return included.reset_index(drop=True), excluded.reset_index(drop=True)


def build(objects: pd.DataFrame, observed: set[tuple[str, str]], max_gap_nt: Optional[int],
          allow_cross_structure_models: bool, allow_exact_duplicates: bool) -> pd.DataFrame:
    records = []
    rows = [r for _, r in objects.iterrows()]
    for idx, (a, b) in enumerate(combinations(rows, 2), start=1):
        g = geometry(int(a.start), int(a.end), int(b.start), int(b.end))
        reasons = []
        if a.species != b.species:
            reasons.append("DIFFERENT_SPECIES")

        both_structure = a.architecture_layer == "RNA_STRUCTURE" and b.architecture_layer == "RNA_STRUCTURE"
        if (both_structure and not allow_cross_structure_models and a.structure_model and b.structure_model
                and a.structure_model != b.structure_model):
            reasons.append("INCOMPATIBLE_STRUCTURE_MODELS")

        duplicate = (
            a.species == b.species and a.pairwise_category == b.pairwise_category
            and int(a.start) == int(b.start) and int(a.end) == int(b.end)
            and a.structure_model == b.structure_model
        )
        if duplicate and not allow_exact_duplicates:
            reasons.append("EXACT_DUPLICATE_ANNOTATION")
        if max_gap_nt is not None and g["relation"] == "SEPARATED" and g["gap_nt"] > max_gap_nt:
            reasons.append("GAP_EXCEEDS_MAXIMUM")

        eligible = not reasons
        homotypy = "HOMOTYPIC" if a.pairwise_category == b.pairwise_category else "HETEROTYPIC"
        layer_relation = "INTRALAYER" if a.architecture_layer == b.architecture_layer else "INTERLAYER"
        pair_cat = " || ".join(sorted((a.pairwise_category, b.pairwise_category)))
        pair_ids = tuple(sorted((a.architecture_id, b.architecture_id)))

        rec = {
            "pair_id": f"BPU_{idx:09d}",
            "pair_key": " || ".join(pair_ids),
            "object_a": a.architecture_id,
            "object_b": b.architecture_id,
            "eligible": eligible,
            "exclusion_reason": ";".join(reasons),
            "eligibility_type": "INELIGIBLE" if not eligible else f"{homotypy}_{layer_relation}_{g['relation']}",
            "homotypy": homotypy,
            "layer_relation": layer_relation,
            "pair_category": pair_cat,
            "observed_interface": pair_ids in observed,
            "a_species": a.species,
            "b_species": b.species,
            "a_category": a.pairwise_category,
            "b_category": b.pairwise_category,
            "a_layer": a.architecture_layer,
            "b_layer": b.architecture_layer,
            "a_structure_model": a.structure_model,
            "b_structure_model": b.structure_model,
            "a_start": int(a.start), "a_end": int(a.end), "a_length_nt": int(a.length_nt),
            "b_start": int(b.start), "b_end": int(b.end), "b_length_nt": int(b.length_nt),
            "center_distance_nt": float(abs(a.center - b.center)),
            "length_difference_nt": int(abs(a.length_nt - b.length_nt)),
        }
        rec.update(g)
        records.append(rec)
    return pd.DataFrame(records)


def main(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw = read_table(args.architecture)
    objects, excluded_objects = prepare(raw)
    observed = load_observed_pairs(args.observed_pairs)
    bpu = build(objects, observed, args.max_gap_nt, args.allow_cross_structure_models,
                args.allow_exact_duplicate_annotations)

    expected = len(objects) * (len(objects) - 1) // 2
    if len(bpu) != expected or bpu["pair_key"].duplicated().any():
        raise RuntimeError("BPU completeness/uniqueness check failed.")

    eligible = bpu[bpu.eligible].copy()
    ineligible = bpu[~bpu.eligible].copy()

    cat = (bpu.groupby(["pair_category", "homotypy", "layer_relation"], dropna=False)
           .agg(total_pairs=("pair_id", "size"), eligible_pairs=("eligible", "sum"),
                observed_interfaces=("observed_interface", "sum"),
                median_gap_nt=("gap_nt", "median"), median_span_nt=("span_length_nt", "median"))
           .reset_index())
    cat["eligible_fraction"] = cat.eligible_pairs / cat.total_pairs
    exclusion = (ineligible.groupby("exclusion_reason", dropna=False).size()
                 .reset_index(name="n_pairs").sort_values("n_pairs", ascending=False))
    geom = (bpu.groupby(["eligible", "relation", "containment"], dropna=False)
            .agg(n_pairs=("pair_id", "size"), median_gap_nt=("gap_nt", "median"),
                 median_boundary_distance_nt=("boundary_distance_nt", "median"),
                 median_span_nt=("span_length_nt", "median"))
            .reset_index())
    obj_cat = (objects.groupby(["species", "pairwise_category", "architecture_layer", "structure_model"], dropna=False)
               .size().reset_index(name="n_objects"))

    outputs = {
        "BPU_all_pairs.csv": bpu,
        "BPU_eligible_pairs.csv": eligible,
        "BPU_ineligible_pairs.csv": ineligible,
        "BPU_object_inventory.csv": objects,
        "BPU_excluded_objects.csv": excluded_objects,
        "BPU_pair_category_summary.csv": cat,
        "BPU_exclusion_summary.csv": exclusion,
        "BPU_geometry_summary.csv": geom,
        "BPU_object_category_summary.csv": obj_cat,
    }
    for name, table in outputs.items():
        table.to_csv(args.output_dir / name, index=False, encoding="utf-8-sig")

    xlsx = args.output_dir / "Boundary_Pair_Universe_29A_1_06A_v1_1.xlsx"
    with pd.ExcelWriter(xlsx, engine="openpyxl") as w:
        if args.full_xlsx:
            bpu.to_excel(w, sheet_name="All_pairs", index=False)
            eligible.to_excel(w, sheet_name="Eligible_pairs", index=False)
            ineligible.to_excel(w, sheet_name="Ineligible_pairs", index=False)
        objects.to_excel(w, sheet_name="Object_inventory", index=False)
        excluded_objects.to_excel(w, sheet_name="Excluded_objects", index=False)
        cat.to_excel(w, sheet_name="Pair_category_summary", index=False)
        exclusion.to_excel(w, sheet_name="Exclusion_summary", index=False)
        geom.to_excel(w, sheet_name="Geometry_summary", index=False)
        obj_cat.to_excel(w, sheet_name="Object_category_summary", index=False)

    qc = {
        "stage": STAGE, "version": VERSION,
        "architecture": str(args.architecture), "architecture_sha256": sha256(args.architecture),
        "n_source_objects": int(len(raw)), "n_pairwise_enabled_qc_pass": int(len(objects)),
        "n_excluded_objects": int(len(excluded_objects)),
        "expected_unordered_pairs": int(expected), "n_pairs_built": int(len(bpu)),
        "n_eligible": int(eligible.shape[0]), "n_ineligible": int(ineligible.shape[0]),
        "eligible_fraction": float(bpu.eligible.mean()),
        "n_pair_categories": int(bpu.pair_category.nunique()),
        "n_observed_interfaces": int(bpu.observed_interface.sum()),
        "n_observed_ineligible": int((bpu.observed_interface & ~bpu.eligible).sum()),
        "max_gap_nt": args.max_gap_nt,
        "allow_cross_structure_models": args.allow_cross_structure_models,
        "allow_exact_duplicate_annotations": args.allow_exact_duplicate_annotations,
        "status": "PASS",
    }
    (args.output_dir / "BPU_QC.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    lines = ["=" * 90, f"{STAGE} v{VERSION} — PASS", "=" * 90] + [f"{k}: {v}" for k, v in qc.items()]
    (args.output_dir / "BPU_QC.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"Workbook: {xlsx}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="29A.1.06A v1.1 — Boundary Pair Universe")
    p.add_argument("--architecture", required=True, type=Path)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--observed-pairs", type=Path, default=None)
    p.add_argument("--max-gap-nt", type=int, default=None)
    p.add_argument("--allow-cross-structure-models", action="store_true")
    p.add_argument("--allow-exact-duplicate-annotations", action="store_true")
    p.add_argument("--full-xlsx", action="store_true")
    a = p.parse_args()
    if not a.architecture.exists(): p.error(f"Not found: {a.architecture}")
    if a.observed_pairs and not a.observed_pairs.exists(): p.error(f"Not found: {a.observed_pairs}")
    if a.max_gap_nt is not None and a.max_gap_nt < 0: p.error("--max-gap-nt must be >= 0")
    return a


if __name__ == "__main__":
    main(parse_args())
