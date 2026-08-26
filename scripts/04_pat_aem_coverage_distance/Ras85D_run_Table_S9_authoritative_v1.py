#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Authoritative Table S9 runner for Ras85D 3′UTR.

Reproduces the final independent logic of:
- Coverage_Enrichment_Test_v2_coverage_only.py
- Functional_Distance_Test_v1.py

Input: authoritative FINAL D. melanogaster rows from 29A_3_01_final_object_registry.csv.
Coordinates: 1-based inclusive; explicit 3′UTR bounds 1..964.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

SEED = 20260718
N_PERM = 10000
UTR_START = 1
UTR_END = 964
APA_WINDOW = 25
SPECIES = "D.melanogaster"

STRUCTURE_CLASSES = ["HELIX", "HAIRPIN", "INTERIOR_LOOP", "MULTILOOP", "BULGE", "EXTERIOR"]
STRUCTURE_MAP = {
    "helices": "HELIX", "helix": "HELIX", "stems": "HELIX", "stem": "HELIX",
    "interior_loops": "INTERIOR_LOOP", "interior_loop": "INTERIOR_LOOP",
    "internal_loops": "INTERIOR_LOOP", "internal_loop": "INTERIOR_LOOP",
    "multiloops": "MULTILOOP", "multiloop": "MULTILOOP",
    "multi_loops": "MULTILOOP", "multi_loop": "MULTILOOP",
    "hairpin_loops": "HAIRPIN", "hairpin_loop": "HAIRPIN",
    "hairpins": "HAIRPIN", "hairpin": "HAIRPIN",
    "bulges": "BULGE", "bulge": "BULGE",
    "exterior": "EXTERIOR", "external": "EXTERIOR", "exterior_loops": "EXTERIOR",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def norm(value: object) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip().lower().replace("′", "'").replace("’", "'")
    return re.sub(r"[\s\-]+", "_", text)


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.fillna("").astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    clean = sorted((int(a), int(b)) for a, b in intervals)
    if not clean:
        return []
    merged: List[List[int]] = [[clean[0][0], clean[0][1]]]
    for start, end in clean[1:]:
        if start <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def mask_from_intervals(intervals: Iterable[Tuple[int, int]], lo: int, hi: int) -> np.ndarray:
    mask = np.zeros(hi - lo + 1, dtype=bool)
    for start, end in merge_intervals(intervals):
        start, end = max(start, lo), min(end, hi)
        if start <= end:
            mask[start - lo:end - lo + 1] = True
    return mask


def bh_adjust(pvalues: Sequence[float]) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    result = np.full(len(p), np.nan)
    valid = np.isfinite(p)
    if not valid.any():
        return result
    values = p[valid]
    order = np.argsort(values)
    ranked = values[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    restored = np.empty(n)
    restored[order] = adjusted
    result[np.where(valid)[0]] = restored
    return result


def randomize_intervals(lengths: Sequence[int], lo: int, hi: int, rng: np.random.Generator) -> List[Tuple[int, int]]:
    utr_length = hi - lo + 1
    randomized = []
    for length in lengths:
        length = int(length)
        if length <= 0:
            continue
        if length >= utr_length:
            randomized.append((lo, hi))
        else:
            first = int(rng.integers(lo, hi - length + 2))
            randomized.append((first, first + length - 1))
    return randomized


def prepare_architecture(registry_path: Path, output_dir: Path) -> pd.DataFrame:
    source = pd.read_csv(registry_path, low_memory=False)
    required = {
        "registry_object_id", "species", "object_type", "feature_group", "feature_class",
        "feature_subclass", "start", "end", "authoritative", "harmonization_status",
        "structure_model", "annotation_value", "site_conservation", "ecm_id",
    }
    missing = sorted(required - set(source.columns))
    if missing:
        raise ValueError("Registry missing: " + ", ".join(missing))

    source["authoritative"] = as_bool(source["authoritative"])
    source["harmonization_status"] = source["harmonization_status"].fillna("").astype(str).str.upper()
    final = source[
        source["authoritative"]
        & source["harmonization_status"].eq("FINAL")
        & source["species"].eq(SPECIES)
    ].copy()

    final["start"] = pd.to_numeric(final["start"], errors="raise").astype(int)
    final["end"] = pd.to_numeric(final["end"], errors="raise").astype(int)
    outside = (final["start"] < UTR_START) | (final["end"] > UTR_END) | (final["start"] > final["end"])
    if outside.any():
        raise ValueError("Invalid/out-of-bounds intervals:\n" + final.loc[outside, ["registry_object_id", "start", "end"]].to_string(index=False))

    # Build the canonical schema expected by the final independent scripts.
    architecture = pd.DataFrame({
        "architecture_id": final["registry_object_id"].astype(str),
        "species": final["species"].astype(str),
        "feature_group": final["feature_group"].astype(str),
        "feature_class": final["feature_class"].astype(str),
        "feature_subclass": final["feature_subclass"].fillna("").astype(str),
        "start": final["start"],
        "end": final["end"],
        "structure_model": final["structure_model"].fillna("").astype(str),
        "object_type": final["object_type"].astype(str),
        "annotation_value": final["annotation_value"].fillna("").astype(str),
        "site_conservation": final["site_conservation"].fillna("").astype(str),
        "ecm_id": final["ecm_id"].fillna("").astype(str),
    })

    # Avoid the substring trap: POORLY_CONSERVED_SITE must not be classified as conserved.
    labels = []
    for row in final.itertuples(index=False):
        if row.object_type == "MIRNA_SITE":
            label = "CONSERVED" if str(row.site_conservation).upper() == "CONSERVED_SITE" else "WEAK"
        elif row.object_type == "STRUCTURE_SEGMENT":
            label = f"{row.structure_model} | {row.feature_subclass}"
        elif row.object_type == "ECM":
            label = str(row.ecm_id)
        else:
            label = f"{row.feature_class} | {row.annotation_value}"
        labels.append(label)
    architecture["architecture_label"] = labels

    if architecture["architecture_id"].duplicated().any():
        raise ValueError("Duplicate architecture_id")

    expected = {
        "MIRNA_SITE": 90, "ECM": 17, "UGUA": 11, "PAS": 9,
        "CLEAVAGE": 12, "DSE": 7, "STRUCTURE_SEGMENT": 364,
    }
    observed = architecture["object_type"].value_counts().to_dict()
    mismatch = {key: (expected[key], int(observed.get(key, 0))) for key in expected if int(observed.get(key, 0)) != expected[key]}
    if mismatch:
        raise ValueError(f"Authoritative Dmel inventory mismatch: {mismatch}")

    out_path = output_dir / "Ras85D_architecture_authoritative_Dmel_v1.csv"
    architecture.to_csv(out_path, index=False, encoding="utf-8-sig")
    return architecture


def classify(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["_structure_model"] = out["structure_model"].map(norm).map({"mfe": "MFE", "centroid": "CENTROID"})
    out["_structure_class"] = None
    structural = out["feature_group"].map(norm).eq("rna_structure")
    for index in out.index[structural]:
        subclass = norm(out.at[index, "feature_subclass"])
        joined = "|".join(norm(out.at[index, column]) for column in ["feature_subclass", "feature_class", "architecture_label"])
        canonical = STRUCTURE_MAP.get(subclass)
        if canonical is None:
            for token, value in STRUCTURE_MAP.items():
                if token in joined:
                    canonical = value
                    break
        out.at[index, "_structure_class"] = canonical
    out["_is_ecm"] = (
        out["feature_group"].map(norm).eq("conservation")
        | out["feature_class"].map(norm).eq("ecm")
        | out["feature_subclass"].map(norm).eq("ecm")
    )
    return out


def feature_masks(df: pd.DataFrame) -> Dict[str, pd.Series]:
    group = df["feature_group"].map(norm)
    feature_class = df["feature_class"].map(norm)
    label = df["architecture_label"].map(norm)
    result: Dict[str, pd.Series] = {}
    for name in ["PAS", "CLEAVAGE", "UGUA", "DSE"]:
        result[name] = group.eq("polya_regulation") & feature_class.str.contains(name.lower(), regex=False, na=False)
    mirna = group.eq("mirna_target") | feature_class.eq("mirna_site")
    conserved = label.eq("conserved")
    result["MIRNA_ALL"] = mirna
    result["MIRNA_CONSERVED"] = mirna & conserved
    result["MIRNA_WEAK"] = mirna & ~conserved
    return result


def structural_masks(df: pd.DataFrame, model: str) -> Dict[str, np.ndarray]:
    subset = df[df["_structure_model"].eq(model)]
    result = {}
    for structural_class in STRUCTURE_CLASSES:
        rows = subset[subset["_structure_class"].eq(structural_class)]
        if not rows.empty:
            result[structural_class] = mask_from_intervals(zip(rows["start"], rows["end"]), UTR_START, UTR_END)
    return result


def strict_consensus_masks(mfe: Dict[str, np.ndarray], centroid: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    result = {}
    for structural_class in STRUCTURE_CLASSES:
        if structural_class in mfe and structural_class in centroid:
            consensus = mfe[structural_class] & centroid[structural_class]
            if consensus.any():
                result[structural_class] = consensus
    return result


def clean_pair(first: np.ndarray, second: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    overlap = first & second
    return first & ~overlap, second & ~overlap


def build_regions(df: pd.DataFrame) -> Tuple[Dict[str, np.ndarray], List[Dict[str, str]], pd.DataFrame]:
    regions: Dict[str, np.ndarray] = {}
    comparisons_list = []
    metadata_rows = []
    mfe = structural_masks(df, "MFE")
    centroid = structural_masks(df, "CENTROID")
    model_maps = {"MFE": mfe, "CENTROID": centroid, "CONSENSUS": strict_consensus_masks(mfe, centroid)}

    for model_name, masks in model_maps.items():
        available = [name for name in STRUCTURE_CLASSES if name in masks and masks[name].any()]
        for class_a, class_b in combinations(available, 2):
            mask_a, mask_b = clean_pair(masks[class_a], masks[class_b])
            if not mask_a.any() or not mask_b.any():
                continue
            key_a = f"STRUCTURE|{model_name}|{class_a}|VS|{class_b}"
            key_b = f"STRUCTURE|{model_name}|{class_b}|VS|{class_a}"
            regions[key_a], regions[key_b] = mask_a, mask_b
            comparisons_list.append({
                "domain": "RNA_STRUCTURE", "test_family": f"STRUCTURE_PAIRWISE_{model_name}",
                "structure_model": model_name, "tested_region_key": key_a, "control_region_key": key_b,
                "tested_region": class_a, "control_region": class_b,
            })

    ecm_rows = df[df["_is_ecm"]].copy()
    all_ecm = mask_from_intervals(zip(ecm_rows["start"], ecm_rows["end"]), UTR_START, UTR_END)
    regions["ECM|ALL"], regions["ECM|NON"] = all_ecm, ~all_ecm
    comparisons_list.append({
        "domain": "ECM", "test_family": "ECM_ALL_VS_NON", "structure_model": "",
        "tested_region_key": "ECM|ALL", "control_region_key": "ECM|NON",
        "tested_region": "ECM_ALL", "control_region": "NON_ECM",
    })
    individual = []
    for row in ecm_rows.itertuples():
        label = str(row.architecture_id)
        key = f"ECM|INDIVIDUAL|{label}"
        regions[key] = mask_from_intervals([(int(row.start), int(row.end))], UTR_START, UTR_END)
        individual.append((key, label))
        comparisons_list.append({
            "domain": "ECM", "test_family": "ECM_INDIVIDUAL_VS_NON", "structure_model": "",
            "tested_region_key": key, "control_region_key": "ECM|NON",
            "tested_region": label, "control_region": "NON_ECM",
        })
    for (key_a, label_a), (key_b, label_b) in combinations(individual, 2):
        mask_a, mask_b = clean_pair(regions[key_a], regions[key_b])
        if not mask_a.any() or not mask_b.any():
            continue
        pair_a, pair_b = f"{key_a}|VS|{label_b}", f"{key_b}|VS|{label_a}"
        regions[pair_a], regions[pair_b] = mask_a, mask_b
        comparisons_list.append({
            "domain": "ECM", "test_family": "ECM_PAIRWISE", "structure_model": "",
            "tested_region_key": pair_a, "control_region_key": pair_b,
            "tested_region": label_a, "control_region": label_b,
        })

    fm = feature_masks(df)
    apa_masks = {}
    for landmark in ["PAS", "CLEAVAGE"]:
        rows = df[fm[landmark]]
        intervals = [
            (max(UTR_START, int(row.start) - APA_WINDOW), min(UTR_END, int(row.end) + APA_WINDOW))
            for row in rows.itertuples()
        ]
        apa_masks[landmark] = mask_from_intervals(intervals, UTR_START, UTR_END)
    apa_masks["APA_COMBINED"] = apa_masks["PAS"] | apa_masks["CLEAVAGE"]
    for landmark, proximal in apa_masks.items():
        distal = ~proximal
        pkey, dkey = f"APA|{landmark}|PROXIMAL", f"APA|{landmark}|DISTAL"
        regions[pkey], regions[dkey] = proximal, distal
        comparisons_list.append({
            "domain": "APA_PROXIMITY", "test_family": f"APA_WINDOW_{landmark}_{APA_WINDOW}NT",
            "structure_model": "", "tested_region_key": pkey, "control_region_key": dkey,
            "tested_region": f"{landmark}_PROXIMAL_{APA_WINDOW}NT", "control_region": f"{landmark}_DISTAL",
        })

    for key, mask in regions.items():
        metadata_rows.append({"species": SPECIES, "region_key": key, "region_length": int(mask.sum())})
    return regions, comparisons_list, pd.DataFrame(metadata_rows)


def stable_log_rr(a: np.ndarray, n_test: int, c: np.ndarray, n_control: int) -> np.ndarray:
    return np.log(((a.astype(float) + 0.5) / (n_test + 1.0)) / ((c.astype(float) + 0.5) / (n_control + 1.0)))


def safe_rr(a: int, n_test: int, c: int, n_control: int) -> float:
    if n_test <= 0 or n_control <= 0:
        return np.nan
    p_test, p_control = a / n_test, c / n_control
    if p_control == 0:
        return np.inf if p_test > 0 else np.nan
    return p_test / p_control


def odds_ratio_ci(a: int, b: int, c: int, d: int) -> Tuple[float, float, float]:
    cells = np.array([a, b, c, d], dtype=float)
    if np.any(cells == 0):
        cells += 0.5
    aa, bb, cc, dd = cells
    odds = aa * dd / (bb * cc)
    se = math.sqrt(1 / aa + 1 / bb + 1 / cc + 1 / dd)
    return odds, math.exp(math.log(odds) - 1.96 * se), math.exp(math.log(odds) + 1.96 * se)


def unique_region_counts(permutation_masks: np.ndarray, regions: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """Calculate counts once per distinct nucleotide mask."""
    packed_perm = np.packbits(permutation_masks, axis=1)
    lookup = np.arange(256, dtype=np.uint8)
    pop = np.unpackbits(lookup[:, None], axis=1).sum(axis=1).astype(np.uint8)
    signature_to_counts = {}
    output = {}
    for key, region in regions.items():
        signature = region.tobytes()
        if signature not in signature_to_counts:
            packed_region = np.packbits(region)
            counts = pop[np.bitwise_and(packed_perm, packed_region)].sum(axis=1).astype(np.int32)
            signature_to_counts[signature] = counts
        output[key] = signature_to_counts[signature]
    return output


def coverage_analysis(df: pd.DataFrame, output_dir: Path) -> Tuple[pd.DataFrame, dict, Dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    regions, comparisons_list, region_metadata = build_regions(df)
    comparison_rows = [{
        "species": SPECIES, **comparison,
        "tested_length": int(regions[comparison["tested_region_key"]].sum()),
        "control_length": int(regions[comparison["control_region_key"]].sum()),
    } for comparison in comparisons_list]

    result_rows = []
    selected_nulls = {}
    for feature_name, row_mask in feature_masks(df).items():
        selected = df[row_mask]
        if selected.empty:
            continue
        merged = merge_intervals(zip(selected["start"], selected["end"]))
        observed_mask = mask_from_intervals(merged, UTR_START, UTR_END)
        lengths = [b - a + 1 for a, b in merged]
        permutation_masks = np.zeros((N_PERM, UTR_END - UTR_START + 1), dtype=bool)
        for index in range(N_PERM):
            permutation_masks[index] = mask_from_intervals(
                randomize_intervals(lengths, UTR_START, UTR_END, rng), UTR_START, UTR_END
            )
        counts = unique_region_counts(permutation_masks, regions)

        for comparison in comparisons_list:
            tested = regions[comparison["tested_region_key"]]
            control = regions[comparison["control_region_key"]]
            n_tested, n_control = int(tested.sum()), int(control.sum())
            a = int((observed_mask & tested).sum())
            c = int((observed_mask & control).sum())
            b, d = n_tested - a, n_control - c
            rr = safe_rr(a, n_tested, c, n_control)
            odds, ci_low, ci_high = odds_ratio_ci(a, b, c, d)
            _, fisher_p = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            perm_tested = counts[comparison["tested_region_key"]]
            perm_control = counts[comparison["control_region_key"]]
            observed_stat = float(stable_log_rr(np.array([a]), n_tested, np.array([c]), n_control)[0])
            null_stats = stable_log_rr(perm_tested, n_tested, perm_control, n_control)
            p_two = (1 + int(np.sum(np.abs(null_stats) >= abs(observed_stat)))) / (N_PERM + 1)
            p_enrich = (1 + int(np.sum(null_stats >= observed_stat))) / (N_PERM + 1)
            p_deplete = (1 + int(np.sum(null_stats <= observed_stat))) / (N_PERM + 1)
            null_rr = np.exp(null_stats)
            direction = ("ENRICHED" if (np.isfinite(rr) and rr > 1) else ("DEPLETED" if (np.isfinite(rr) and rr < 1) else ("ENRICHED" if (not np.isfinite(rr) and a > 0) else "UNDEFINED")))
            result_rows.append({
                "species": SPECIES, "feature": feature_name, "domain": comparison["domain"],
                "test_family": comparison["test_family"], "structure_model": comparison["structure_model"],
                "tested_region": comparison["tested_region"], "control_region": comparison["control_region"],
                "n_feature_intervals_merged": len(merged), "tested_length": n_tested, "control_length": n_control,
                "covered_test": a, "covered_control": c, "free_test": b, "free_control": d,
                "coverage_test": a / n_tested if n_tested else np.nan,
                "coverage_control": c / n_control if n_control else np.nan,
                "RR": rr, "OR": odds, "OR_CI95_low": ci_low, "OR_CI95_high": ci_high,
                "fisher_p_two_sided": fisher_p, "observed_log_RR_stat": observed_stat,
                "null_RR_median": float(np.median(null_rr)), "null_RR_q025": float(np.quantile(null_rr, 0.025)),
                "null_RR_q975": float(np.quantile(null_rr, 0.975)), "perm_p_two_sided": p_two,
                "perm_p_enrichment": p_enrich, "perm_p_depletion": p_deplete, "effect_direction": direction,
            })

    results = pd.DataFrame(result_rows)
    results["fdr_family"] = results["species"].astype(str) + "|" + results["test_family"].astype(str) + "|" + results["feature"].astype(str)
    results["FDR_BH"] = np.nan
    for _, indices in results.groupby("fdr_family", dropna=False).groups.items():
        results.loc[indices, "FDR_BH"] = bh_adjust(results.loc[indices, "perm_p_two_sided"].to_numpy())
    results["result"] = "NOT_SIGNIFICANT"
    results.loc[(results["FDR_BH"] < 0.05) & (results["RR"] > 1), "result"] = "SIGNIFICANT_ENRICHMENT"
    results.loc[(results["FDR_BH"] < 0.05) & (results["RR"] < 1), "result"] = "SIGNIFICANT_DEPLETION"
    results = results.sort_values(["species", "domain", "test_family", "feature", "FDR_BH", "perm_p_two_sided"], na_position="last").reset_index(drop=True)

    summary = results.groupby(["species", "domain", "test_family", "feature"], dropna=False).agg(
        n_tests=("feature", "size"), n_nominal=("perm_p_two_sided", lambda x: int((x < 0.05).sum())),
        n_fdr=("FDR_BH", lambda x: int((x < 0.05).sum())), min_perm_p=("perm_p_two_sided", "min"),
        min_FDR=("FDR_BH", "min"), median_RR=("RR", "median")
    ).reset_index()

    results.to_csv(output_dir / "coverage_enrichment_results_v2.csv", index=False, encoding="utf-8-sig")
    results[results.FDR_BH < 0.05].to_csv(output_dir / "coverage_enrichment_significant_v2.csv", index=False, encoding="utf-8-sig")
    results[results.perm_p_two_sided < 0.05].to_csv(output_dir / "coverage_enrichment_nominal_v2.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(output_dir / "coverage_enrichment_family_summary_v2.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(comparison_rows).drop_duplicates().to_csv(output_dir / "coverage_region_comparisons_v2.csv", index=False, encoding="utf-8-sig")
    region_metadata.to_csv(output_dir / "coverage_region_definitions_v2.csv", index=False, encoding="utf-8-sig")

    classified = df.copy()
    for name, mask in feature_masks(df).items():
        classified[f"_feature_{name}"] = mask
    classified.to_csv(output_dir / "architecture_classified_coverage_v2.csv", index=False, encoding="utf-8-sig")

    qc = {
        "test": "Coverage Enrichment Test v2.0 authoritative",
        "n_architecture_rows": int(len(df)), "n_permutations": N_PERM, "apa_window_nt": APA_WINDOW,
        "seed": SEED, "region_comparisons": int(len(comparison_rows)), "tests": int(len(results)),
        "nominal_p_lt_0_05": int((results.perm_p_two_sided < 0.05).sum()),
        "fdr_lt_0_05": int((results.FDR_BH < 0.05).sum()),
        "feature_counts": {name: int(mask.sum()) for name, mask in feature_masks(df).items()},
        "structure_counts": {str(k): int(v) for k, v in df["_structure_class"].value_counts(dropna=False).items()},
        "structure_model_counts": {str(k): int(v) for k, v in df["_structure_model"].value_counts(dropna=False).items()},
        "ecm_rows": int(df["_is_ecm"].sum()),
    }
    (output_dir / "Coverage_Enrichment_Test_v2_QC.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    return results, qc, selected_nulls


def interval_distance(feature: Tuple[int, int], landmark: Tuple[int, int]) -> int:
    fs, fe = feature; ls, le = landmark
    if fe < ls: return ls - fe - 1
    if le < fs: return fs - le - 1
    return 0


def nearest_distances(features: Sequence[Tuple[int, int]], landmarks: Sequence[Tuple[int, int]]) -> np.ndarray:
    return np.array([min(interval_distance(f, l) for l in landmarks) for f in features], dtype=float)


def distance_analysis(df: pd.DataFrame, output_dir: Path) -> Tuple[pd.DataFrame, dict, Dict[str, np.ndarray]]:
    rng = np.random.default_rng(SEED)
    fm = feature_masks(df)
    catalog = {
        "PAS": merge_intervals(zip(df[fm["PAS"]]["start"], df[fm["PAS"]]["end"])),
        "CLEAVAGE": merge_intervals(zip(df[fm["CLEAVAGE"]]["start"], df[fm["CLEAVAGE"]]["end"])),
    }
    landmark_rows = [{
        "species": SPECIES, "landmark": name, "n_intervals_merged": len(intervals),
        "total_landmark_length": int(sum(b - a + 1 for a, b in intervals)),
    } for name, intervals in catalog.items()]

    rows = []
    nulls = {}
    for feature_name, mask in fm.items():
        selected = df[mask]
        if selected.empty:
            continue
        merged = merge_intervals(zip(selected["start"], selected["end"]))
        lengths = [b - a + 1 for a, b in merged]
        perm_sets = [randomize_intervals(lengths, UTR_START, UTR_END, rng) for _ in range(N_PERM)]
        for landmark_name, landmarks in catalog.items():
            if feature_name == landmark_name:
                continue
            observed = nearest_distances(merged, landmarks)
            null_medians = np.array([np.median(nearest_distances(intervals, landmarks)) for intervals in perm_sets], dtype=float)
            null_means = np.array([np.mean(nearest_distances(intervals, landmarks)) for intervals in perm_sets], dtype=float)
            null_zeros = np.array([np.mean(nearest_distances(intervals, landmarks) == 0) for intervals in perm_sets], dtype=float)
            observed_median = float(np.median(observed))
            null_center = float(np.median(null_medians))
            deviation = observed_median - null_center
            p_two = (1 + int(np.sum(np.abs(null_medians - null_center) >= abs(deviation)))) / (N_PERM + 1)
            p_closer = (1 + int(np.sum(null_medians <= observed_median))) / (N_PERM + 1)
            p_farther = (1 + int(np.sum(null_medians >= observed_median))) / (N_PERM + 1)
            direction = "CLOSER_TO_LANDMARK" if observed_median < null_center else ("FARTHER_FROM_LANDMARK" if observed_median > null_center else "NEUTRAL")
            key = f"{feature_name}__{landmark_name}"
            nulls[key] = null_medians
            rows.append({
                "species": SPECIES, "feature": feature_name, "landmark": landmark_name,
                "n_feature_intervals_merged": len(merged), "n_landmark_intervals_merged": len(landmarks),
                "observed_median_distance_nt": observed_median, "observed_mean_distance_nt": float(np.mean(observed)),
                "observed_overlap_fraction": float(np.mean(observed == 0)), "null_median_distance_nt": null_center,
                "null_median_q025": float(np.quantile(null_medians, 0.025)), "null_median_q975": float(np.quantile(null_medians, 0.975)),
                "null_mean_distance_nt": float(np.median(null_means)), "null_overlap_fraction": float(np.median(null_zeros)),
                "median_distance_shift_nt": deviation, "perm_p_two_sided": p_two, "perm_p_closer": p_closer,
                "perm_p_farther": p_farther, "effect_direction": direction,
            })

    results = pd.DataFrame(rows)
    results["fdr_family"] = results["species"].astype(str) + "|DISTANCE_PANEL|" + results["feature"].astype(str)
    results["FDR_BH"] = np.nan
    for _, indices in results.groupby("fdr_family", dropna=False).groups.items():
        results.loc[indices, "FDR_BH"] = bh_adjust(results.loc[indices, "perm_p_two_sided"].to_numpy())
    results["result"] = "NOT_SIGNIFICANT"
    results.loc[(results.FDR_BH < 0.05) & results.effect_direction.eq("CLOSER_TO_LANDMARK"), "result"] = "SIGNIFICANTLY_CLOSER"
    results.loc[(results.FDR_BH < 0.05) & results.effect_direction.eq("FARTHER_FROM_LANDMARK"), "result"] = "SIGNIFICANTLY_FARTHER"
    results = results.sort_values(["species", "feature", "FDR_BH", "perm_p_two_sided", "landmark"]).reset_index(drop=True)
    summary = results.groupby(["species", "feature"], dropna=False).agg(
        n_landmarks_tested=("landmark", "size"), n_nominal=("perm_p_two_sided", lambda x: int((x < 0.05).sum())),
        n_fdr=("FDR_BH", lambda x: int((x < 0.05).sum())), min_perm_p=("perm_p_two_sided", "min"), min_FDR=("FDR_BH", "min")
    ).reset_index()

    results.to_csv(output_dir / "functional_distance_results_v1.csv", index=False, encoding="utf-8-sig")
    results[results.FDR_BH < 0.05].to_csv(output_dir / "functional_distance_significant_v1.csv", index=False, encoding="utf-8-sig")
    results[results.perm_p_two_sided < 0.05].to_csv(output_dir / "functional_distance_nominal_v1.csv", index=False, encoding="utf-8-sig")
    summary.to_csv(output_dir / "functional_distance_summary_v1.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(landmark_rows).to_csv(output_dir / "functional_distance_landmarks_v1.csv", index=False, encoding="utf-8-sig")
    np.savez_compressed(output_dir / "functional_distance_null_medians_v1.npz", **nulls)

    classified = df.copy()
    for name, mask in fm.items():
        classified[f"_feature_{name}"] = mask
    classified.to_csv(output_dir / "architecture_classified_distance_v1.csv", index=False, encoding="utf-8-sig")
    qc = {
        "test": "Functional Distance Test v1.0 authoritative", "n_architecture_rows": int(len(df)),
        "n_permutations": N_PERM, "seed": SEED, "landmarks_argument": "PAS,CLEAVAGE", "exclude_self": True,
        "tests": int(len(results)), "nominal_p_lt_0_05": int((results.perm_p_two_sided < 0.05).sum()),
        "fdr_lt_0_05": int((results.FDR_BH < 0.05).sum()),
    }
    (output_dir / "Functional_Distance_Test_v1_QC.json").write_text(json.dumps(qc, indent=2, ensure_ascii=False), encoding="utf-8")
    return results, qc, nulls


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.output_root
    if root.exists():
        import shutil
        shutil.rmtree(root)
    coverage_dir = root / "Coverage_Enrichment_Test_v2_authoritative"
    distance_dir = root / "Functional_Distance_Test_v1_authoritative"
    coverage_dir.mkdir(parents=True)
    distance_dir.mkdir(parents=True)

    architecture = prepare_architecture(args.registry, root)
    classified = classify(architecture)

    coverage_results, coverage_qc, _ = coverage_analysis(classified, coverage_dir)
    distance_results, distance_qc, _ = distance_analysis(classified, distance_dir)

    summary = {
        "status": "PASS", "registry": str(args.registry), "registry_sha256": sha256(args.registry),
        "architecture_rows": int(len(architecture)), "coverage": coverage_qc, "distance": distance_qc,
        "coverage_significant_rows": coverage_results[coverage_results.FDR_BH < 0.05].to_dict(orient="records"),
        "distance_significant_rows": distance_results[distance_results.FDR_BH < 0.05].to_dict(orient="records"),
    }
    (root / "Ras85D_Table_S9_authoritative_run_QC_v1.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({
        "status": "PASS", "architecture_rows": len(architecture),
        "coverage_tests": len(coverage_results), "coverage_nominal": coverage_qc["nominal_p_lt_0_05"],
        "coverage_fdr": coverage_qc["fdr_lt_0_05"], "distance_tests": len(distance_results),
        "distance_nominal": distance_qc["nominal_p_lt_0_05"], "distance_fdr": distance_qc["fdr_lt_0_05"],
    }, indent=2))


if __name__ == "__main__":
    main()
