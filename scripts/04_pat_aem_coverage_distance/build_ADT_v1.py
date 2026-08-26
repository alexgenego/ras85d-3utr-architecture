#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Build Ras85D Architecture Domain Table v1 (ADT v1).

Inputs (auto-discovered under /content):
  Ras85D_Master_Architecture_Table_v1_1.csv
  Ras85D_IEL_master_table_final.csv (or Ras85D_IEL_master_table.csv)
  arpip_34_regions_alignment3_species_specific.csv

Outputs in /content:
  Ras85D_Architecture_Domain_Table_v1_long.csv       (IEL × species; 102 rows)
  Ras85D_Architecture_Domain_Table_v1_consensus.csv  (one row per IEL; 34 rows)
  Ras85D_Architecture_Domain_Table_v1.xlsx
  Ras85D_ADT_v1_QC_report.txt

Important conventions:
- ECM_BLOCK intervals are the exact curated ECM mask.
- ECM_OCCURRENCE rows are MAST anchors, not exact motif boundaries.
- Cleavage sites are boundaries between two alignment columns.
- MFE and centroid structural models are summarized separately.
"""

from __future__ import annotations

import glob
import os
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

BASE = Path('/content')
SPECIES_PREFIX = {
    'D.melanogaster': 'dmel',
    'D.yakuba': 'dyak',
    'D.virilis': 'dvir',
}
SPECIES = list(SPECIES_PREFIX)
APA_NEAR = 25
APA_CORE = 15
ECM_INTERFACE = 10
STRUCT_BOUNDARY = 2


def newest(patterns: Sequence[str]) -> Optional[Path]:
    hits: List[str] = []
    for p in patterns:
        hits.extend(glob.glob(p, recursive=True))
    hits = sorted(set(hits), key=os.path.getmtime, reverse=True)
    return Path(hits[0]) if hits else None


def required(label: str, patterns: Sequence[str]) -> Path:
    path = newest(patterns)
    if path is None:
        raise FileNotFoundError(f'Не найден файл: {label}\n' + '\n'.join(patterns))
    return path


MAT_PATH = required('MAT v1.1', [
    '/content/Ras85D_Master_Architecture_Table_v1_1.csv',
    '/content/**/*Ras85D_Master_Architecture_Table_v1_1.csv',
])
IEL_PATH = required('IEL master', [
    '/content/Ras85D_IEL_master_table_final.csv',
    '/content/Ras85D_IEL_master_table.csv',
    '/content/**/*Ras85D_IEL_master_table_final.csv',
    '/content/**/*Ras85D_IEL_master_table.csv',
])
PROJ_PATH = required('IEL projections', [
    '/content/arpip_34_regions_alignment3_species_specific.csv',
    '/content/**/*arpip_34_regions_alignment3_species_specific.csv',
])
HOTSPOT_PATH = newest([
    '/content/**/*Ras85D_spatial_deletion_hotspot_analysis.xlsx',
    '/content/drive/MyDrive/**/*Ras85D_spatial_deletion_hotspot_analysis.xlsx',
])


def txt(v: Any) -> str:
    return '' if pd.isna(v) else str(v).strip()


def intval(v: Any) -> Optional[int]:
    try:
        return None if pd.isna(v) else int(float(v))
    except (TypeError, ValueError):
        return None


def join_unique(values: Iterable[Any]) -> str:
    return ';'.join(sorted({txt(v) for v in values if txt(v) not in {'', 'nan', '<NA>'}}))


def getv(row: pd.Series, names: Sequence[str], default: Any = pd.NA) -> Any:
    for name in names:
        if name in row.index and pd.notna(row[name]):
            return row[name]
    return default


def overlap_len(a: int, b: int, c: int, d: int) -> int:
    return max(0, min(b, d) - max(a, c) + 1)


def interval_dist(a: int, b: int, c: int, d: int) -> int:
    if b < c:
        return c - b
    if d < a:
        return a - d
    return 0


def relation(a: int, b: int, c: int, d: int) -> str:
    ov = overlap_len(a, b, c, d)
    if ov == 0:
        return 'upstream' if b < c else 'downstream'
    if a == c and b == d:
        return 'coincident'
    if c <= a and b <= d:
        if a == c:
            return 'inside_at_left_boundary'
        if b == d:
            return 'inside_at_right_boundary'
        return 'inside'
    if a <= c and d <= b:
        return 'contains'
    if a < c <= b:
        return 'crosses_left_boundary'
    if a <= d < b:
        return 'crosses_right_boundary'
    return 'partial_overlap'


def structure_model(source: Any) -> str:
    s = txt(source).lower()
    if 'centroid' in s:
        return 'centroid'
    if 'mfe' in s:
        return 'MFE'
    return 'unspecified'


def structure_type(value: Any) -> str:
    s = txt(value).lower().replace('-', '_').replace(' ', '_')
    mapping = {
        'helices': 'helix', 'helix': 'helix',
        'hairpin_loops': 'hairpin_loop', 'hairpin': 'hairpin_loop',
        'interior_loops': 'internal_loop', 'internal_loops': 'internal_loop',
        'interior_loop': 'internal_loop', 'internal_loop': 'internal_loop',
        'multiloops': 'multiloop', 'multiloop': 'multiloop',
        'bulges': 'bulge', 'bulge': 'bulge',
        'external_loops': 'external_loop', 'external_loop': 'external_loop',
    }
    return mapping.get(s, s or 'unspecified')


mat = pd.read_csv(MAT_PATH)
iel = pd.read_csv(IEL_PATH)
proj = pd.read_csv(PROJ_PATH)

if 'IEL_ID' not in iel.columns:
    iel['IEL_ID'] = [f'IEL_{i+1:03d}' for i in range(len(iel))]
if 'region_id' not in iel.columns:
    iel['region_id'] = np.arange(len(iel))
if 'region_id' not in proj.columns:
    raise ValueError('В таблице проекций нет region_id.')

# Add hotspot summaries when available.
if HOTSPOT_PATH:
    try:
        xls = pd.ExcelFile(HOTSPOT_PATH)
        for sheet in ['deletion_region_mrca_summary', 'region_geometry_summary']:
            if sheet not in xls.sheet_names:
                continue
            extra = pd.read_excel(HOTSPOT_PATH, sheet_name=sheet)
            if extra['region_id'].min() == 1 and iel['region_id'].min() == 0 and extra['region_id'].max() == len(iel):
                extra['region_id'] -= 1
            iel = iel.merge(extra, on='region_id', how='left', suffixes=('', f'_{sheet}'))
    except Exception as exc:
        print('Предупреждение hotspot Excel:', exc)

full = iel.merge(proj, on='region_id', how='left', suffixes=('', '_projection'))

# MAT layers.
def layer(name: str, species: Optional[str] = None) -> pd.DataFrame:
    out = mat[mat['object_type'] == name].copy()
    return out if species is None else out[out['species'] == species].copy()

ECM = layer('ECM_BLOCK')
MAST = layer('ECM_OCCURRENCE')
UGUA = layer('UGUA')
PAS = layer('PAS')
CLEAVAGE = layer('CLEAVAGE_SITE')
DSE = layer('DOWNSTREAM_ELEMENT')
STRUCT = layer('RNA_STRUCTURE')
STRUCT['structure_model'] = STRUCT['source_column'].apply(structure_model)
STRUCT['structure_type_norm'] = STRUCT['object_subtype'].apply(structure_type)


def projection_interval(row: pd.Series, prefix: str) -> Dict[str, Any]:
    starts = [f'{prefix}_alignment3_start', f'{prefix}_start_alignment3', f'{prefix}_start']
    ends = [f'{prefix}_alignment3_end', f'{prefix}_end_alignment3', f'{prefix}_end']
    labels = [f'{prefix}_alignment3_label', f'{prefix}_location_label', f'{prefix}_alignment3_location_label']
    classes = [f'{prefix}_projection_class', 'projection_class']
    s, e = intval(getv(row, starts)), intval(getv(row, ends))
    cls = txt(getv(row, classes, 'unknown')) or 'unknown'
    if s is not None and e is not None:
        return {'start': min(s,e), 'end': max(s,e), 'label': f'{min(s,e)}-{max(s,e)}', 'support': 'direct_interval', 'class': cls}
    label = next((txt(row[n]) for n in labels if n in row.index and txt(row[n])), '')
    nums = [int(x) for x in re.findall(r'\d+', label)]
    if '|' in label and len(nums) >= 2:
        return {'start': min(nums[:2]), 'end': max(nums[:2]), 'label': label, 'support': 'between_columns', 'class': cls}
    if '-' in label and len(nums) >= 2:
        return {'start': min(nums[:2]), 'end': max(nums[:2]), 'label': label, 'support': 'label_interval', 'class': cls}
    if len(nums) == 1:
        return {'start': nums[0], 'end': nums[0], 'label': label, 'support': 'label_point', 'class': cls}
    return {'start': None, 'end': None, 'label': label, 'support': 'unresolved', 'class': cls}


def summarize_intervals(qs: int, qe: int, objects: pd.DataFrame, prefix: str) -> Dict[str, Any]:
    if objects.empty:
        return {f'{prefix}_n_overlap': 0, f'{prefix}_nearest_distance': pd.NA}
    w = objects.copy()
    w['s'] = pd.to_numeric(w['alignment3_start'], errors='coerce')
    w['e'] = pd.to_numeric(w['alignment3_end'], errors='coerce')
    w = w.dropna(subset=['s','e']).copy()
    w[['s','e']] = w[['s','e']].astype(int)
    w['overlap_length'] = w.apply(lambda r: overlap_len(qs,qe,r.s,r.e), axis=1)
    w['distance'] = w.apply(lambda r: interval_dist(qs,qe,r.s,r.e), axis=1)
    w['relation'] = w.apply(lambda r: relation(qs,qe,r.s,r.e), axis=1)
    ov = w[w.overlap_length > 0]
    near = w.sort_values(['distance','s','e','object_id']).iloc[0]
    left = w[w.e < qs].sort_values(['e','s'])
    right = w[w.s > qe].sort_values(['s','e'])
    out = {
        f'{prefix}_n_overlap': len(ov),
        f'{prefix}_overlap_ids': join_unique(ov.object_id),
        f'{prefix}_overlap_subtypes': join_unique(ov.object_subtype),
        f'{prefix}_overlap_relations': join_unique(ov.relation),
        f'{prefix}_overlap_total_length': int(ov.overlap_length.sum()) if len(ov) else 0,
        f'{prefix}_nearest_id': near.object_id,
        f'{prefix}_nearest_subtype': near.object_subtype,
        f'{prefix}_nearest_distance': int(near.distance),
        f'{prefix}_nearest_relation': near.relation,
    }
    if len(left):
        x = left.iloc[-1]
        out.update({f'{prefix}_left_id': x.object_id, f'{prefix}_left_subtype': x.object_subtype, f'{prefix}_left_distance': qs-int(x.e)})
    else:
        out.update({f'{prefix}_left_id': pd.NA, f'{prefix}_left_subtype': pd.NA, f'{prefix}_left_distance': pd.NA})
    if len(right):
        x = right.iloc[0]
        out.update({f'{prefix}_right_id': x.object_id, f'{prefix}_right_subtype': x.object_subtype, f'{prefix}_right_distance': int(x.s)-qe})
    else:
        out.update({f'{prefix}_right_id': pd.NA, f'{prefix}_right_subtype': pd.NA, f'{prefix}_right_distance': pd.NA})
    return out


def summarize_cleavage(qs: int, qe: int, objects: pd.DataFrame) -> Dict[str, Any]:
    if objects.empty:
        return {'cleavage_n_overlap': 0, 'cleavage_nearest_distance': pd.NA}
    w = objects.copy()
    w['left'] = pd.to_numeric(w['alignment3_boundary_left'], errors='coerce').fillna(pd.to_numeric(w['alignment3_start'], errors='coerce'))
    w['right'] = pd.to_numeric(w['alignment3_boundary_right'], errors='coerce').fillna(pd.to_numeric(w['alignment3_end'], errors='coerce'))
    w = w.dropna(subset=['left','right']).copy()
    w[['left','right']] = w[['left','right']].astype(int)
    w['overlap'] = w.apply(lambda r: not (qe < r.left or qs > r.right), axis=1)
    w['distance'] = w.apply(lambda r: 0 if r.overlap else min(abs(qs-r.left),abs(qs-r.right),abs(qe-r.left),abs(qe-r.right)), axis=1)
    ov = w[w.overlap]
    near = w.sort_values(['distance','left','right']).iloc[0]
    center = (near.left+near.right)/2
    qcenter = (qs+qe)/2
    side = 'overlap' if near.distance == 0 else ('upstream' if qcenter < center else 'downstream')
    return {
        'cleavage_n_overlap': len(ov),
        'cleavage_overlap_ids': join_unique(ov.object_id),
        'cleavage_overlap_subtypes': join_unique(ov.object_subtype),
        'cleavage_nearest_id': near.object_id,
        'cleavage_nearest_subtype': near.object_subtype,
        'cleavage_nearest_distance': int(near.distance),
        'cleavage_nearest_side': side,
        'cleavage_nearest_boundary': f'{int(near.left)}|{int(near.right)}',
    }


def summarize_mast(qs: int, qe: int, objects: pd.DataFrame) -> Dict[str, Any]:
    if objects.empty:
        return {'mast_nearest_distance': pd.NA, 'mast_n_within_25': 0}
    w = objects.copy()
    w['p'] = pd.to_numeric(w['alignment3_start'], errors='coerce')
    w = w.dropna(subset=['p']).copy(); w['p'] = w.p.astype(int)
    w['distance'] = w.p.apply(lambda p: 0 if qs <= p <= qe else min(abs(p-qs), abs(p-qe)))
    near = w.sort_values(['distance','p']).iloc[0]
    around = w[w.distance <= 25]
    return {
        'mast_nearest_id': near.object_id,
        'mast_nearest_family': near.motif_family,
        'mast_nearest_distance': int(near.distance),
        'mast_nearest_orientation': near.orientation,
        'mast_n_within_25': len(around),
        'mast_families_within_25': join_unique(around.motif_family),
        'mast_orientations_within_25': join_unique(around.orientation),
        'mast_reverse_within_25': bool((around.orientation == 'reverse_complement').any()),
    }


def summarize_structure(qs: int, qe: int, objects: pd.DataFrame, model: str) -> Dict[str, Any]:
    key = f'structure_{model.lower()}'
    w = objects[objects.structure_model == model].copy()
    if w.empty:
        return {f'{key}_n_overlap': 0, f'{key}_primary_type': pd.NA, f'{key}_nearest_distance': pd.NA}
    w['s'] = pd.to_numeric(w.alignment3_start, errors='coerce')
    w['e'] = pd.to_numeric(w.alignment3_end, errors='coerce')
    w = w.dropna(subset=['s','e']).copy(); w[['s','e']] = w[['s','e']].astype(int)
    w['overlap_length'] = w.apply(lambda r: overlap_len(qs,qe,r.s,r.e), axis=1)
    w['distance'] = w.apply(lambda r: interval_dist(qs,qe,r.s,r.e), axis=1)
    w['relation'] = w.apply(lambda r: relation(qs,qe,r.s,r.e), axis=1)
    ov = w[w.overlap_length > 0]
    near = w.sort_values(['distance','s','e']).iloc[0]
    if len(ov):
        primary = ov.sort_values(['overlap_length','s'], ascending=[False,True]).iloc[0]
        ptype, prel = primary.structure_type_norm, primary.relation
        bd = min(min(abs(qs-r.s),abs(qs-r.e),abs(qe-r.s),abs(qe-r.e)) for _,r in ov.iterrows())
    else:
        ptype, prel = pd.NA, pd.NA
        bd = min(abs(qs-near.s),abs(qs-near.e),abs(qe-near.s),abs(qe-near.e))
    return {
        f'{key}_n_overlap': len(ov),
        f'{key}_overlap_types': join_unique(ov.structure_type_norm),
        f'{key}_overlap_ids': join_unique(ov.object_id),
        f'{key}_primary_type': ptype,
        f'{key}_primary_relation': prel,
        f'{key}_nearest_type': near.structure_type_norm,
        f'{key}_nearest_id': near.object_id,
        f'{key}_nearest_distance': int(near.distance),
        f'{key}_boundary_distance': int(bd),
        f'{key}_near_boundary': bool(bd <= STRUCT_BOUNDARY),
    }


def ecm_context(r: Dict[str, Any]) -> str:
    if r.get('ecm_n_overlap',0) > 0:
        rel = txt(r.get('ecm_overlap_relations'))
        if 'coincident' in rel: return 'coincident_with_ECM'
        if 'crosses' in rel: return 'crosses_ECM_boundary'
        if 'inside' in rel: return 'inside_ECM'
        if 'contains' in rel: return 'IEL_contains_ECM'
        return 'overlaps_ECM'
    ld, rd = intval(r.get('ecm_left_distance')), intval(r.get('ecm_right_distance'))
    if ld is not None and rd is not None:
        return 'ECM_interface' if min(ld,rd) <= ECM_INTERFACE else 'between_ECM'
    if ld is None and rd is not None: return 'before_first_ECM'
    if ld is not None and rd is None: return 'after_last_ECM'
    return 'ECM_unresolved'


def apa_context(r: Dict[str, Any]) -> Tuple[str,str]:
    feats = []
    for key,label in [('ugua_nearest_distance','UGUA'),('pas_nearest_distance','PAS'),('cleavage_nearest_distance','CLEAVAGE'),('dse_nearest_distance','DSE')]:
        d = intval(r.get(key))
        if d is not None and d <= APA_NEAR: feats.append(label)
    s = set(feats)
    if {'UGUA','PAS','CLEAVAGE','DSE'} <= s: cls = 'complete_APA_neighborhood'
    elif {'PAS','CLEAVAGE','DSE'} <= s: cls = 'core_APA_neighborhood'
    elif {'PAS','CLEAVAGE'} <= s: cls = 'PAS_cleavage_neighborhood'
    elif {'CLEAVAGE','DSE'} <= s: cls = 'cleavage_DSE_neighborhood'
    elif feats: cls = 'partial_APA_neighborhood'
    else: cls = 'APA_distant'
    return '+'.join(feats) if feats else 'none', cls


def structure_consensus(r: Dict[str, Any]) -> Tuple[str,str]:
    m, c = txt(r.get('structure_mfe_primary_type')), txt(r.get('structure_centroid_primary_type'))
    if m and c:
        return (m, 'MFE_centroid_agreement') if m == c else (f'MFE:{m}|centroid:{c}', 'MFE_centroid_disagreement')
    if m: return m, 'MFE_only'
    if c: return c, 'centroid_only'
    return 'none', 'no_structure_overlap'


def domain_class(r: Dict[str, Any]) -> str:
    ec = r['ECM_context']; ap = r['APA_context_class']; st = r['structure_consensus_type']
    has_ecm = ec in {'inside_ECM','coincident_with_ECM','crosses_ECM_boundary','overlaps_ECM','IEL_contains_ECM'}
    interface = ec in {'ECM_interface','between_ECM'}
    has_apa = ap != 'APA_distant'
    has_struct = st not in {'','none'}
    cd, pd, dd = intval(r.get('cleavage_nearest_distance')), intval(r.get('pas_nearest_distance')), intval(r.get('dse_nearest_distance'))
    apa_core = cd is not None and cd <= APA_CORE and ((pd is not None and pd <= APA_CORE) or (dd is not None and dd <= APA_CORE))
    if has_ecm and has_apa and has_struct: return 'composite_ECM_APA_structure_domain'
    if interface and has_apa and has_struct: return 'composite_interface_APA_structure_domain'
    if has_ecm and has_apa: return 'ECM_APA_domain'
    if has_ecm and has_struct: return 'ECM_structure_domain'
    if interface and has_apa: return 'ECM_interface_APA_domain'
    if interface and has_struct: return 'ECM_interface_structure_domain'
    if apa_core and has_struct: return 'structural_APA_core'
    if apa_core: return 'APA_core'
    if has_apa and has_struct: return 'APA_structure_domain'
    if has_apa: return 'APA_neighborhood'
    if has_ecm: return 'ECM_domain'
    if interface: return 'ECM_interface'
    if has_struct: return 'structure_only_domain'
    return 'unclassified_architectural_context'


rows: List[Dict[str,Any]] = []
for _, ir in full.iterrows():
    for species in SPECIES:
        pfx = SPECIES_PREFIX[species]
        pr = projection_interval(ir,pfx)
        r: Dict[str,Any] = {
            'IEL_ID': ir.IEL_ID,
            'region_id': intval(ir.region_id),
            'species': species,
            'alignment37_start': getv(ir,['start_alignment','alignment37_start']),
            'alignment37_end': getv(ir,['end_alignment','alignment37_end']),
            'alignment37_length': getv(ir,['length','region_length','alignment37_length']),
            'IEL_alignment3_start': pr['start'], 'IEL_alignment3_end': pr['end'],
            'IEL_alignment3_length': pr['end']-pr['start']+1 if pr['start'] is not None else pd.NA,
            'IEL_alignment3_label': pr['label'], 'IEL_projection_class': pr['class'],
            'IEL_support_type': pr['support'], 'IEL_projection_resolved': pr['start'] is not None,
            'n_deletion_events': getv(ir,['n_deletion_events','n_overlapping_deletion_events','n_local_deletion_events']),
            'n_insertion_events': getv(ir,['n_insertion_events','n_overlapping_insertion_events','n_local_insertion_events']),
            'n_all_pairs': getv(ir,['n_all_pairs']),
            'n_informative_topological_pairs': getv(ir,['n_informative_geometry_pairs','n_informative_topological_pairs']),
            'n_independent_recurrent_pairs': getv(ir,['n_independent_recurrent_pairs']),
            'n_probable_ancestral_remodeling_pairs': getv(ir,['n_probable_ancestral_remodeling_pairs']),
            'n_ambiguous_ancestral_pairs': getv(ir,['n_ambiguous_ancestral_pairs']),
            'independent_fraction': getv(ir,['independent_fraction']),
            'remodeling_fraction': getv(ir,['remodeling_fraction']),
            'IEL_local_composition': getv(ir,['IEL_local_composition','indel_composition']),
            'IEL_topological_hotspot_class': getv(ir,['hotspot_geometry_class','hotspot_topological_class','IEL_topological_class']),
        }
        if pr['start'] is None:
            r.update({'ECM_context':'unresolved','APA_feature_profile':'unresolved','APA_context_class':'unresolved','structure_consensus_type':'unresolved','structure_model_agreement':'unresolved','architectural_formula':'unresolved','architectural_domain_class':'unresolved'})
            rows.append(r); continue
        qs,qe = pr['start'],pr['end']
        r.update(summarize_intervals(qs,qe,ECM[ECM.species==species],'ecm'))
        r.update(summarize_mast(qs,qe,MAST[MAST.species==species]))
        r.update(summarize_intervals(qs,qe,UGUA[UGUA.species==species],'ugua'))
        r.update(summarize_intervals(qs,qe,PAS[PAS.species==species],'pas'))
        r.update(summarize_cleavage(qs,qe,CLEAVAGE[CLEAVAGE.species==species]))
        r.update(summarize_intervals(qs,qe,DSE[DSE.species==species],'dse'))
        r.update(summarize_structure(qs,qe,STRUCT[STRUCT.species==species],'MFE'))
        r.update(summarize_structure(qs,qe,STRUCT[STRUCT.species==species],'centroid'))
        r['ECM_context'] = ecm_context(r)
        r['APA_feature_profile'], r['APA_context_class'] = apa_context(r)
        r['structure_consensus_type'], r['structure_model_agreement'] = structure_consensus(r)
        comps=[]
        if r['ECM_context'] in {'inside_ECM','coincident_with_ECM','crosses_ECM_boundary','overlaps_ECM','IEL_contains_ECM'}: comps.append('ECM')
        elif r['ECM_context'] in {'ECM_interface','between_ECM'}: comps.append('ECM_INTERFACE')
        if r['APA_context_class'] != 'APA_distant': comps.append('APA')
        if r['structure_consensus_type'] != 'none': comps.append('RNA_STRUCTURE')
        r['architectural_formula'] = '+'.join(comps) if comps else 'UNANNOTATED'
        r['architectural_domain_class'] = domain_class(r)
        rows.append(r)

adt_long = pd.DataFrame(rows)

# Consensus table: one row per IEL.
consensus=[]
for iid,g in adt_long.groupby('IEL_ID',sort=False):
    f=g.iloc[0]; resolved=g[g.IEL_projection_resolved==True]
    classes=sorted(set(resolved.architectural_domain_class.dropna().astype(str)))
    if not classes: cross='unresolved'
    elif len(classes)==1: cross='conserved_domain_class'
    elif len(classes)==2: cross='partially_conserved_domain_class'
    else: cross='species_variable_domain_class'
    out={
        'IEL_ID':iid,'region_id':f.region_id,'alignment37_start':f.alignment37_start,'alignment37_end':f.alignment37_end,'alignment37_length':f.alignment37_length,
        'n_deletion_events':f.get('n_deletion_events',pd.NA),'n_insertion_events':f.get('n_insertion_events',pd.NA),'n_independent_recurrent_pairs':f.get('n_independent_recurrent_pairs',pd.NA),'independent_fraction':f.get('independent_fraction',pd.NA),
        'IEL_local_composition':f.get('IEL_local_composition',pd.NA),'IEL_topological_hotspot_class':f.get('IEL_topological_hotspot_class',pd.NA),
        'n_species_resolved':len(resolved),'species_resolved':join_unique(resolved.species),'cross_species_domain_consensus':cross,
        'architectural_domain_classes':';'.join(classes),'architectural_formulas':join_unique(resolved.architectural_formula),
        'ECM_contexts':join_unique(resolved.ECM_context),'APA_context_classes':join_unique(resolved.APA_context_class),'structure_consensus_types':join_unique(resolved.structure_consensus_type),
        'n_species_ECM_overlap':int((resolved.get('ecm_n_overlap',0).fillna(0).astype(float)>0).sum()) if len(resolved) else 0,
        'n_species_APA_near':int((resolved.APA_context_class!='APA_distant').sum()),
        'n_species_structure_overlap':int((resolved.structure_consensus_type!='none').sum()),
    }
    for species,short in [('D.melanogaster','Dmel'),('D.yakuba','Dyak'),('D.virilis','Dvir')]:
        sub=g[g.species==species]
        if len(sub):
            x=sub.iloc[0]
            out[f'{short}_alignment3_label']=x.IEL_alignment3_label
            out[f'{short}_ECM_context']=x.ECM_context
            out[f'{short}_APA_context']=x.APA_context_class
            out[f'{short}_structure_type']=x.structure_consensus_type
            out[f'{short}_domain_class']=x.architectural_domain_class
    consensus.append(out)
adt_consensus=pd.DataFrame(consensus).sort_values('region_id').reset_index(drop=True)

OUT_LONG=BASE/'Ras85D_Architecture_Domain_Table_v1_long.csv'
OUT_CONS=BASE/'Ras85D_Architecture_Domain_Table_v1_consensus.csv'
OUT_XLSX=BASE/'Ras85D_Architecture_Domain_Table_v1.xlsx'
OUT_QC=BASE/'Ras85D_ADT_v1_QC_report.txt'
adt_long.to_csv(OUT_LONG,index=False)
adt_consensus.to_csv(OUT_CONS,index=False)

with pd.ExcelWriter(OUT_XLSX,engine='openpyxl') as w:
    adt_consensus.to_excel(w,sheet_name='ADT_consensus_34_IEL',index=False)
    adt_long.to_excel(w,sheet_name='ADT_species_102_rows',index=False)
    adt_long.architectural_domain_class.value_counts(dropna=False).rename_axis('class').reset_index(name='n').to_excel(w,sheet_name='domain_class_summary',index=False)
    adt_long.ECM_context.value_counts(dropna=False).rename_axis('ECM_context').reset_index(name='n').to_excel(w,sheet_name='ECM_context_summary',index=False)
    adt_long.APA_context_class.value_counts(dropna=False).rename_axis('APA_context').reset_index(name='n').to_excel(w,sheet_name='APA_context_summary',index=False)
    adt_long.structure_consensus_type.value_counts(dropna=False).rename_axis('structure').reset_index(name='n').to_excel(w,sheet_name='structure_summary',index=False)

qc='\n'.join([
    'Ras85D ADT v1 QC REPORT','='*60,
    f'MAT: {MAT_PATH}',f'IEL: {IEL_PATH}',f'Projection: {PROJ_PATH}',
    f'ADT long rows: {len(adt_long)} (expected {len(iel)*3})',
    f'ADT consensus rows: {len(adt_consensus)} (expected {len(iel)})','',
    'Resolved projections:',adt_long.groupby('species').IEL_projection_resolved.agg(['sum','count']).to_string(),'',
    'Domain classes:',adt_long.architectural_domain_class.value_counts(dropna=False).to_string(),'',
    'ECM contexts:',adt_long.ECM_context.value_counts(dropna=False).to_string(),'',
    'APA contexts:',adt_long.APA_context_class.value_counts(dropna=False).to_string(),'',
    'Structure agreement:',adt_long.structure_model_agreement.value_counts(dropna=False).to_string(),
])
OUT_QC.write_text(qc,encoding='utf-8')

print('\nADT v1 СОЗДАНА')
print('Long:',OUT_LONG,adt_long.shape)
print('Consensus:',OUT_CONS,adt_consensus.shape)
print('Excel:',OUT_XLSX)
print('QC:',OUT_QC)
print('\nАрхитектурные классы:')
print(adt_long.architectural_domain_class.value_counts(dropna=False))
print('\nECM-контекст:')
print(adt_long.ECM_context.value_counts(dropna=False))
print('\nAPA-контекст:')
print(adt_long.APA_context_class.value_counts(dropna=False))
print('\nСогласие MFE/centroid:')
print(adt_long.structure_model_agreement.value_counts(dropna=False))
