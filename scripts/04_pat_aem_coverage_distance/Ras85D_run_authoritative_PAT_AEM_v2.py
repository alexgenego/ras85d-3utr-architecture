#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run the original PAT v2 and build_ADT_v1 scripts on the final manuscript-associated registry.

Required files:
- 29A_3_01_final_object_registry.csv
- Ras85D_IEL_master_table_final.csv
- arpip_34_regions_alignment3_species_specific.csv
- Pairwise_Architecture_Test_v2.py
- build_ADT_v1.py

The wrapper patches only file paths. Statistical code remains the original code.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(command: list[str], cwd: Path) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def patch_pat(source: Path, destination: Path, workspace: Path, pat_out: Path) -> None:
    text = source.read_text(encoding="utf-8")
    substitutions = {
        r"BACKUP_ROOT\s*=\s*Path\([^\n]+\)":
            f"BACKUP_ROOT = Path({str(workspace)!r})",
        r"OUTPUT_DIR\s*=\s*Path\([^\n]+\)":
            f"OUTPUT_DIR = Path({str(pat_out)!r})",
        r"BACKUP_DIR\s*=\s*BACKUP_ROOT\s*/\s*'Pairwise_Architecture_Test_v2'":
            f"BACKUP_DIR = Path({str(pat_out / 'Backup')!r})",
        r"BAN_V2_BACKUP_DIR\s*=\s*BACKUP_ROOT\s*/\s*'Basal_Architecture_Network_v2'":
            f"BAN_V2_BACKUP_DIR = Path({str(pat_out / 'BAN_v2')!r})",
    }
    for pattern, replacement in substitutions.items():
        text, count = re.subn(pattern, replacement, text, count=1)
        if count != 1:
            raise RuntimeError(f"PAT patch failed for: {pattern}")
    destination.write_text(text, encoding="utf-8")


def patch_adt(source: Path, destination: Path, workspace: Path) -> None:
    text = source.read_text(encoding="utf-8")
    text, count = re.subn(
        r"BASE\s*=\s*Path\([^\n]+\)",
        f"BASE = Path({str(workspace)!r})",
        text,
        count=1,
    )
    if count != 1:
        raise RuntimeError("ADT patch failed for BASE")

    # The historical ADT script auto-discovers files using literal /content paths.
    # Redirect only those filesystem paths to the prepared workspace; statistical
    # calculations and classification logic remain unchanged.
    text = text.replace("'/content", f"'{str(workspace)}")
    text = text.replace('"/content', f'"{str(workspace)}')
    destination.write_text(text, encoding="utf-8")


def validate_outputs(root: Path) -> dict:
    pat_dir = root / "PAT_authoritative"
    aem_dir = root / "AEM_authoritative"

    pat_inventory_path = pat_dir / "Ras85D_Pairwise_Architecture_Test_v2_inventory.csv"
    pat_results_path = pat_dir / "Ras85D_Pairwise_Architecture_Test_v2_results.csv"
    pat_edges_path = pat_dir / "Ras85D_Pairwise_Architecture_Test_v2_edge_summary.csv"

    required = [pat_inventory_path, pat_results_path, pat_edges_path]
    for path in required:
        if not path.exists():
            raise FileNotFoundError(f"Missing PAT output: {path}")

    inventory = pd.read_csv(pat_inventory_path)
    results = pd.read_csv(pat_results_path)
    edges = pd.read_csv(pat_edges_path)

    expected = {
        "D.melanogaster": (17, 39, 11, 9, 12, 7, 364),
        "D.yakuba": (17, 36, 9, 9, 8, 10, 334),
        "D.virilis": (15, 52, 12, 13, 9, 18, 437),
    }
    for species, values in expected.items():
        row = inventory[inventory["species"].eq(species)]
        if len(row) != 1:
            raise ValueError(f"Missing/duplicate PAT inventory row for {species}")
        row = row.iloc[0]
        observed = (
            int(row["n_ECM"]),
            int(row["n_APA"]),
            int(row["n_UGUA"]),
            int(row["n_PAS"]),
            int(row["n_CLEAVAGE"]),
            int(row["n_DSE"]),
            int(row["n_MFE_structure"] + row["n_centroid_structure"]),
        )
        if observed != values:
            raise ValueError(
                f"PAT inventory mismatch for {species}: "
                f"expected={values}; observed={observed}"
            )

    if len(results) != 43:
        raise ValueError(f"Expected 43 PAT metrics, found {len(results)}")
    if len(edges) != 3:
        raise ValueError(f"Expected 3 PAT edge rows, found {len(edges)}")

    long_path = aem_dir / "Ras85D_Architecture_Domain_Table_v1_long.csv"
    consensus_path = aem_dir / "Ras85D_Architecture_Domain_Table_v1_consensus.csv"
    if not long_path.exists() or not consensus_path.exists():
        raise FileNotFoundError("Missing rebuilt AEM/ADT outputs")

    long = pd.read_csv(long_path)
    consensus = pd.read_csv(consensus_path)
    if len(long) != 102:
        raise ValueError(f"Expected 102 AEM long rows, found {len(long)}")
    if len(consensus) != 34:
        raise ValueError(f"Expected 34 AEM consensus rows, found {len(consensus)}")

    return {
        "status": "PASS",
        "PAT_metrics": int(len(results)),
        "PAT_edges": int(len(edges)),
        "PAT_inventory": inventory.to_dict(orient="records"),
        "AEM_long_rows": int(len(long)),
        "AEM_consensus_rows": int(len(consensus)),
        "PAT_edge_summary": edges.to_dict(orient="records"),
        "AEM_consensus_class_counts": (
            consensus["cross_species_domain_consensus"]
            .value_counts(dropna=False)
            .to_dict()
            if "cross_species_domain_consensus" in consensus.columns
            else {}
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--iel-master", type=Path, required=True)
    parser.add_argument("--iel-projections", type=Path, required=True)
    parser.add_argument("--pat-script", type=Path, required=True)
    parser.add_argument("--adt-script", type=Path, required=True)
    parser.add_argument("--prepare-script", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()

    for path in [
        args.registry, args.iel_master, args.iel_projections,
        args.pat_script, args.adt_script, args.prepare_script,
    ]:
        if not path.exists():
            raise FileNotFoundError(path)

    root = args.output_root.resolve()
    workspace = root / "Workspace"
    pat_out = root / "PAT_authoritative"
    aem_out = root / "AEM_authoritative"

    if root.exists():
        shutil.rmtree(root)
    workspace.mkdir(parents=True)
    pat_out.mkdir(parents=True)
    aem_out.mkdir(parents=True)

    run(
        [
            sys.executable,
            str(args.prepare_script.resolve()),
            "--registry", str(args.registry.resolve()),
            "--output-dir", str(workspace),
            "--alignment-length", "1315",
        ],
        cwd=workspace,
    )

    shutil.copy2(args.iel_master, workspace / "Ras85D_IEL_master_table_final.csv")
    shutil.copy2(
        args.iel_projections,
        workspace / "arpip_34_regions_alignment3_species_specific.csv",
    )

    patched_pat = workspace / "Pairwise_Architecture_Test_v2_authoritative.py"
    patched_adt = workspace / "build_ADT_v1_authoritative.py"
    patch_pat(args.pat_script, patched_pat, workspace, pat_out)
    patch_adt(args.adt_script, patched_adt, workspace)

    run([sys.executable, str(patched_pat)], cwd=workspace)

    # build_ADT writes into BASE/workspace. Move its final outputs into AEM.
    run([sys.executable, str(patched_adt)], cwd=workspace)
    for name in [
        "Ras85D_Architecture_Domain_Table_v1_long.csv",
        "Ras85D_Architecture_Domain_Table_v1_consensus.csv",
        "Ras85D_Architecture_Domain_Table_v1.xlsx",
        "Ras85D_ADT_v1_QC_report.txt",
    ]:
        source = workspace / name
        if source.exists():
            shutil.move(str(source), aem_out / name)

    summary = validate_outputs(root)
    summary["input_sha256"] = {
        "registry": sha256(args.registry),
        "iel_master": sha256(args.iel_master),
        "iel_projections": sha256(args.iel_projections),
        "pat_script": sha256(args.pat_script),
        "adt_script": sha256(args.adt_script),
    }
    summary["output_root"] = str(root)
    qc_path = root / "Ras85D_PAT_AEM_authoritative_rerun_QC_v1.json"
    qc_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    print("=" * 78)
    print("AUTHORITATIVE PAT/AEM RERUN — PASS")
    print("=" * 78)
    print("PAT metrics:", summary["PAT_metrics"])
    print("PAT edges:", summary["PAT_edges"])
    print("AEM long rows:", summary["AEM_long_rows"])
    print("AEM consensus rows:", summary["AEM_consensus_rows"])
    print("QC:", qc_path)


if __name__ == "__main__":
    main()
