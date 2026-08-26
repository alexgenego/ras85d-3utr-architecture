#!/usr/bin/env python3
"""Rebuild Data S1 coordinate maps from the corrected 2276-column FASTA.

Expected input: Data_S1_02_alignment_names_corrected_2276cols.fasta
Coordinate convention: 1-based inclusive.
Analytical matrix: alignment columns 1-2247.
"""
from pathlib import Path
import csv, gzip

def read_fasta(path):
    rec=[]; name=None; parts=[]
    for raw in Path(path).read_text(encoding="utf-8-sig").splitlines():
        line=raw.strip()
        if not line: continue
        if line.startswith(">"):
            if name is not None: rec.append((name,"".join(parts)))
            name=line[1:].strip(); parts=[]
        else: parts.append(line)
    if name is not None: rec.append((name,"".join(parts)))
    return rec

records=read_fasta("Data_S1_02_alignment_names_corrected_2276cols.fasta")
assert len(records)==37 and {len(s) for _,s in records}=={2276}
with gzip.open("rebuilt_alignment_to_taxon_coordinates_long.tsv.gz","wt",encoding="utf-8",newline="") as out:
    w=csv.writer(out,delimiter="\t",lineterminator="\n")
    w.writerow(["taxon","alignment_column_1based","symbol","is_gap","taxon_coordinate_full_1based","taxon_coordinate_analytical_1based","in_analytical_matrix","in_excluded_distal_segment"])
    for taxon,seq in records:
        cfull=canal=0
        for i,ch in enumerate(seq,1):
            gap=ch=="-"
            if not gap: cfull+=1
            if i<=2247 and not gap: canal+=1
            w.writerow([taxon,i,ch,int(gap),"" if gap else cfull,"" if gap or i>2247 else canal,int(i<=2247),int(i>=2248)])
print("Coordinate map rebuilt.")
