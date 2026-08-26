#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
29A.1.06B — Boundary Pair Universe permutation generator

Consumes BPU_eligible_pairs.csv produced by 29A.1.06A v1.1.

Null model
----------
For every permutation, sample without replacement from the eligible universe
while preserving the observed number of interfaces within matching strata:
species, homotypy, layer_relation, relation, gap_bin, span_bin and
length_pair_bin.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

STAGE = "29A.1.06B"
VERSION = "1.0"
TRUE_VALUES = {"true", "t", "1", "yes", "y"}


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    return series.astype(str).str.strip().str.lower().isin(TRUE_VALUES)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    values = p[ok]
    if not len(values):
        return out
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.minimum(adjusted, 1.0)
    out[np.where(ok)[0]] = restored
    return out


def quantile_bin(series: pd.Series, reference: pd.Series, q: int, label: str) -> pd.Series:
    ref = pd.to_numeric(reference, errors="coerce").dropna()
    if ref.nunique() < 2:
        return pd.Series([f"{label}_ALL"] * len(series), index=series.index)
    edges = np.unique(np.quantile(ref, np.linspace(0, 1, q + 1)))
    if len(edges) < 2:
        return pd.Series([f"{label}_ALL"] * len(series), index=series.index)
    edges[0], edges[-1] = -np.inf, np.inf
    return pd.cut(
        pd.to_numeric(series, errors="coerce"),
        bins=edges,
        include_lowest=True,
        duplicates="drop",
    ).astype(str).radd(f"{label}_")


def attach_observed(bpu: pd.DataFrame, observed_path: Path | None) -> pd.DataFrame:
    result = bpu.copy()
    if observed_path is None:
        if "observed_interface" not in result.columns:
            raise ValueError("No observed_interface column and no observed-pairs table.")
        result["observed_interface"] = as_bool(result["observed_interface"])
        return result

    observed = (
        pd.read_csv(observed_path)
        if observed_path.suffix.lower() == ".csv"
        else pd.read_excel(observed_path)
    )
    alternatives = [
        ("object_a", "object_b"),
        ("architecture_id_a", "architecture_id_b"),
        ("source", "target"),
        ("node_a", "node_b"),
    ]
    columns = next((item for item in alternatives if set(item).issubset(observed.columns)), None)
    if columns is None:
        raise ValueError("Unsupported observed-pairs columns.")

    keys = {
        " || ".join(sorted((str(a).strip(), str(b).strip())))
        for a, b in observed[list(columns)].itertuples(index=False, name=None)
        if str(a).strip() != str(b).strip()
    }
    result["observed_interface"] = result["pair_key"].isin(keys)
    return result


def build_strata(data: pd.DataFrame, n_bins: int) -> tuple[pd.DataFrame, list[str]]:
    observed = data[data.observed_interface]
    data = data.copy()
    data["species"] = np.where(
        data.a_species.eq(data.b_species),
        data.a_species,
        "CROSS_SPECIES",
    )
    data["gap_bin"] = quantile_bin(data.gap_nt, observed.gap_nt, n_bins, "GAP")
    data["span_bin"] = quantile_bin(data.span_length_nt, observed.span_length_nt, n_bins, "SPAN")
    mean_length = (data.a_length_nt + data.b_length_nt) / 2
    observed_mean_length = (observed.a_length_nt + observed.b_length_nt) / 2
    data["length_pair_bin"] = quantile_bin(
        mean_length,
        observed_mean_length,
        n_bins,
        "LEN",
    )
    strata = [
        "species",
        "homotypy",
        "layer_relation",
        "relation",
        "gap_bin",
        "span_bin",
        "length_pair_bin",
    ]
    data["stratum"] = data[strata].astype(str).agg(" | ".join, axis=1)
    return data, strata


def feasibility(data: pd.DataFrame) -> pd.DataFrame:
    return (
        data.groupby("stratum", dropna=False)
        .agg(
            eligible_pairs=("pair_id", "size"),
            observed_interfaces=("observed_interface", "sum"),
        )
        .reset_index()
        .assign(feasible=lambda frame: frame.observed_interfaces <= frame.eligible_pairs)
    )


def summarize_counts(sample: pd.DataFrame, categories: list[str]) -> np.ndarray:
    counts = sample.pair_category.value_counts()
    return np.array([counts.get(category, 0) for category in categories], dtype=np.int32)


def main(args: argparse.Namespace) -> None:
    args.output_dir.mkdir(parents=True, exist_ok=True)
    bpu = pd.read_csv(args.bpu)
    required = {
        "pair_id", "pair_key", "eligible", "pair_category", "homotypy",
        "layer_relation", "relation", "gap_nt", "span_length_nt",
        "a_length_nt", "b_length_nt", "a_species", "b_species",
    }
    missing = sorted(required - set(bpu.columns))
    if missing:
        raise ValueError("BPU missing columns: " + ", ".join(missing))

    bpu["eligible"] = as_bool(bpu["eligible"])
    bpu = bpu[bpu.eligible].copy().reset_index(drop=True)
    bpu = attach_observed(bpu, args.observed_pairs)
    if not bpu.observed_interface.any():
        raise ValueError("Zero observed interfaces.")

    bpu, strata_columns = build_strata(bpu, args.geometry_bins)
    feasible = feasibility(bpu)
    if not feasible.feasible.all():
        raise RuntimeError(
            "Infeasible strata:\n" +
            feasible[~feasible.feasible].to_string(index=False)
        )

    observed = bpu[bpu.observed_interface].copy()
    categories = sorted(bpu.pair_category.unique())
    observed_counts = summarize_counts(observed, categories)

    strata_groups = []
    for stratum, group in bpu.groupby("stratum", sort=False):
        count = int(group.observed_interface.sum())
        if count:
            strata_groups.append((stratum, group.index.to_numpy(), count))

    rng = np.random.default_rng(args.seed)
    null = np.zeros((args.n_permutations, len(categories)), dtype=np.int32)
    global_metrics = [
        "HOMOTYPIC", "HETEROTYPIC", "INTRALAYER", "INTERLAYER",
        "OVERLAP", "TOUCHING", "SEPARATED",
    ]
    global_null = {
        metric: np.zeros(args.n_permutations, dtype=np.int32)
        for metric in global_metrics
    }

    checkpoint_every = max(1, args.checkpoint_every)
    for permutation_index in range(args.n_permutations):
        chosen_parts = []
        for _, indices, count in strata_groups:
            chosen_parts.append(rng.choice(indices, size=count, replace=False))
        chosen = np.concatenate(chosen_parts)
        sample = bpu.loc[chosen]
        null[permutation_index, :] = summarize_counts(sample, categories)

        for metric in ["HOMOTYPIC", "HETEROTYPIC"]:
            global_null[metric][permutation_index] = int(sample.homotypy.eq(metric).sum())
        for metric in ["INTRALAYER", "INTERLAYER"]:
            global_null[metric][permutation_index] = int(sample.layer_relation.eq(metric).sum())
        for metric in ["OVERLAP", "TOUCHING", "SEPARATED"]:
            global_null[metric][permutation_index] = int(sample.relation.eq(metric).sum())

        if (
            (permutation_index + 1) % checkpoint_every == 0
            or permutation_index + 1 == args.n_permutations
        ):
            np.savez_compressed(
                args.output_dir / "BPU_permutation_checkpoint.npz",
                completed=permutation_index + 1,
                categories=np.array(categories, dtype=object),
                null=null[:permutation_index + 1],
            )
            print(
                f"Permutations completed: "
                f"{permutation_index + 1}/{args.n_permutations}",
                flush=True,
            )

    rows = []
    for column, category in enumerate(categories):
        array = null[:, column]
        observed_value = int(observed_counts[column])
        null_mean = float(array.mean())
        p_enrichment = float(
            (1 + np.sum(array >= observed_value)) /
            (args.n_permutations + 1)
        )
        p_depletion = float(
            (1 + np.sum(array <= observed_value)) /
            (args.n_permutations + 1)
        )
        rows.append({
            "pair_category": category,
            "observed": observed_value,
            "null_mean": null_mean,
            "null_sd": float(array.std(ddof=1)),
            "null_q025": float(np.quantile(array, 0.025)),
            "null_median": float(np.median(array)),
            "null_q975": float(np.quantile(array, 0.975)),
            "observed_minus_expected": observed_value - null_mean,
            "fold_enrichment": (
                observed_value / null_mean
                if null_mean > 0
                else (np.inf if observed_value > 0 else np.nan)
            ),
            "p_enrichment": p_enrichment,
            "p_depletion": p_depletion,
            "p_two_sided": min(1.0, 2 * min(p_enrichment, p_depletion)),
        })

    results = pd.DataFrame(rows)
    results["fdr_two_sided"] = bh(results.p_two_sided.to_numpy())
    results["direction"] = np.where(
        results.observed > results.null_mean,
        "ENRICHED",
        np.where(results.observed < results.null_mean, "DEPLETED", "EQUAL"),
    )
    results["result"] = np.where(
        results.fdr_two_sided < 0.05,
        "SIGNIFICANT_" + results.direction,
        "NOT_SIGNIFICANT",
    )
    results = results.sort_values(
        ["fdr_two_sided", "p_two_sided", "pair_category"]
    )

    global_rows = []
    for metric in global_metrics:
        if metric in {"HOMOTYPIC", "HETEROTYPIC"}:
            observed_value = int(observed.homotypy.eq(metric).sum())
        elif metric in {"INTRALAYER", "INTERLAYER"}:
            observed_value = int(observed.layer_relation.eq(metric).sum())
        else:
            observed_value = int(observed.relation.eq(metric).sum())

        array = global_null[metric]
        global_rows.append({
            "metric": metric,
            "observed": observed_value,
            "null_mean": float(array.mean()),
            "null_sd": float(array.std(ddof=1)),
            "null_q025": float(np.quantile(array, 0.025)),
            "null_q975": float(np.quantile(array, 0.975)),
            "p_two_sided": min(
                1.0,
                2 * min(
                    (1 + np.sum(array >= observed_value)) /
                    (args.n_permutations + 1),
                    (1 + np.sum(array <= observed_value)) /
                    (args.n_permutations + 1),
                ),
            ),
        })
    global_results = pd.DataFrame(global_rows)

    stratum_summary = (
        bpu.groupby("stratum", dropna=False)
        .agg(
            eligible_pairs=("pair_id", "size"),
            observed_interfaces=("observed_interface", "sum"),
            n_pair_categories=("pair_category", "nunique"),
        )
        .reset_index()
    )

    observed.to_csv(
        args.output_dir / "BPU_observed_interfaces_used.csv",
        index=False,
        encoding="utf-8-sig",
    )
    results.to_csv(
        args.output_dir / "BPU_permutation_pair_category_results.csv",
        index=False,
        encoding="utf-8-sig",
    )
    global_results.to_csv(
        args.output_dir / "BPU_permutation_global_results.csv",
        index=False,
        encoding="utf-8-sig",
    )
    stratum_summary.to_csv(
        args.output_dir / "BPU_permutation_strata.csv",
        index=False,
        encoding="utf-8-sig",
    )
    np.savez_compressed(
        args.output_dir / "BPU_null_distributions.npz",
        categories=np.array(categories, dtype=object),
        null_counts=null,
        n_permutations=args.n_permutations,
        seed=args.seed,
    )

    with pd.ExcelWriter(
        args.output_dir / "Boundary_Pair_Permutation_29A_1_06B.xlsx",
        engine="openpyxl",
    ) as writer:
        results.to_excel(writer, "Pair_category_results", index=False)
        results[results.fdr_two_sided < 0.05].to_excel(
            writer,
            "FDR_significant",
            index=False,
        )
        global_results.to_excel(writer, "Global_results", index=False)
        stratum_summary.to_excel(writer, "Strata", index=False)
        feasible.to_excel(writer, "Feasibility", index=False)
        observed.to_excel(writer, "Observed_interfaces", index=False)

    qc = {
        "stage": STAGE,
        "version": VERSION,
        "status": "PASS",
        "bpu": str(args.bpu),
        "bpu_sha256": sha256(args.bpu),
        "n_eligible_pairs": int(len(bpu)),
        "n_observed_interfaces": int(len(observed)),
        "n_strata_with_observed": int(len(strata_groups)),
        "n_pair_categories": int(len(categories)),
        "n_permutations": int(args.n_permutations),
        "seed": int(args.seed),
        "geometry_bins": int(args.geometry_bins),
        "strata_columns": strata_columns,
        "n_fdr_significant_pair_categories": int(
            (results.fdr_two_sided < 0.05).sum()
        ),
    }
    (args.output_dir / "BPU_permutation_QC.json").write_text(
        json.dumps(qc, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(qc, indent=2, ensure_ascii=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="29A.1.06B — BPU permutation generator"
    )
    parser.add_argument("--bpu", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--observed-pairs", type=Path, default=None)
    parser.add_argument("--n-permutations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260718)
    parser.add_argument("--geometry-bins", type=int, default=4)
    parser.add_argument("--checkpoint-every", type=int, default=100)
    args = parser.parse_args()
    if not args.bpu.exists():
        parser.error(f"Not found: {args.bpu}")
    if args.observed_pairs and not args.observed_pairs.exists():
        parser.error(f"Not found: {args.observed_pairs}")
    if args.n_permutations < 1:
        parser.error("--n-permutations must be >= 1")
    if args.geometry_bins < 1:
        parser.error("--geometry-bins must be >= 1")
    return args


if __name__ == "__main__":
    main(parse_args())
