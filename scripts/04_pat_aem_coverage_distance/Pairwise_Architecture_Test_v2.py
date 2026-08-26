#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Ras85D Pairwise Architecture Test v2
=====================================

Автономный анализ трёх рёбер базовой архитектурной сети:
1. ECM <-> RNA secondary structure
2. ECM <-> APA architecture
3. RNA secondary structure <-> APA architecture

Нулевая модель: полный перебор всех ненулевых циклических сдвигов
одного слоя относительно второго вдоль alignment3.
"""

from __future__ import annotations

import shutil
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd

BACKUP_ROOT = Path('/content/drive/MyDrive/Ras85D_ADT_backup')
MAT_CANDIDATES = [
    BACKUP_ROOT / 'Ras85D_Master_Architecture_Table_v1_1.csv',
    BACKUP_ROOT / 'Ras85D_Master_Architecture_Table_v1.csv',
]
MASTER3_CANDIDATES = [BACKUP_ROOT / 'Ras85D_alignment3_master_base.csv']
BAN_EDGE_CANDIDATES = [BACKUP_ROOT / 'Basal_Architecture_Network_v1' / 'Ras85D_BAN_v1_edges.csv']
BAN_NODE_CANDIDATES = [BACKUP_ROOT / 'Basal_Architecture_Network_v1' / 'Ras85D_BAN_v1_nodes.csv']
OUTPUT_DIR = Path('/content')
BACKUP_DIR = BACKUP_ROOT / 'Pairwise_Architecture_Test_v2'
BAN_V2_BACKUP_DIR = BACKUP_ROOT / 'Basal_Architecture_Network_v2'
SPECIES_ORDER = ['D.melanogaster', 'D.yakuba', 'D.virilis']
APA_TYPES = ['UGUA', 'PAS', 'CLEAVAGE_SITE', 'DOWNSTREAM_ELEMENT']
BOUNDARY_WINDOW = 2

Segment = Tuple[int, int]


def first_existing(candidates: Sequence[Path], label: str) -> Path:
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f'Не найден {label}. Проверены:\n' + '\n'.join(map(str, candidates)))


def normalize_segment(start: int, end: int) -> Segment:
    return (min(start, end), max(start, end))


def interval_overlap(query_segments: Sequence[Segment], target_intervals: np.ndarray) -> bool:
    if len(target_intervals) == 0:
        return False
    for qs, qe in query_segments:
        overlap = np.minimum(qe, target_intervals[:, 1]) - np.maximum(qs, target_intervals[:, 0]) + 1
        if np.any(overlap > 0):
            return True
    return False


def minimum_interval_distance(query_segments: Sequence[Segment], target_intervals: np.ndarray) -> float:
    if len(target_intervals) == 0:
        return np.inf
    minimum = np.inf
    for qs, qe in query_segments:
        distances = np.where(
            qe < target_intervals[:, 0],
            target_intervals[:, 0] - qe,
            np.where(target_intervals[:, 1] < qs, qs - target_intervals[:, 1], 0),
        )
        minimum = min(minimum, float(distances.min()))
    return minimum


def minimum_boundary_distance(query_segments: Sequence[Segment], boundaries: np.ndarray) -> float:
    if len(boundaries) == 0:
        return np.inf
    minimum = np.inf
    for qs, qe in query_segments:
        inside = (boundaries >= qs) & (boundaries <= qe)
        if inside.any():
            return 0.0
        distances = np.minimum(np.abs(boundaries - qs), np.abs(boundaries - qe))
        minimum = min(minimum, float(distances.min()))
    return minimum


def shift_position_circular(position: int, shift: int, sequence_length: int) -> int:
    return ((int(position) - 1 + int(shift)) % int(sequence_length)) + 1


def shift_segment_circular(start: int, end: int, shift: int, sequence_length: int) -> List[Segment]:
    ss = shift_position_circular(start, shift, sequence_length)
    se = shift_position_circular(end, shift, sequence_length)
    if start == end:
        return [(ss, ss)]
    if ss <= se:
        return [(ss, se)]
    return [(ss, sequence_length), (1, se)]


def shift_segments_circular(segments: Sequence[Segment], shift: int, sequence_length: int) -> List[Segment]:
    shifted: List[Segment] = []
    for start, end in segments:
        shifted.extend(shift_segment_circular(start, end, shift, sequence_length))
    return shifted


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p_values = np.asarray(p_values, dtype=float)
    n = len(p_values)
    if n == 0:
        return np.asarray([], dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted_ranked = ranked * n / np.arange(1, n + 1)
    adjusted_ranked = np.minimum.accumulate(adjusted_ranked[::-1])[::-1]
    adjusted_ranked = np.minimum(adjusted_ranked, 1.0)
    adjusted = np.empty(n, dtype=float)
    adjusted[order] = adjusted_ranked
    return adjusted


def add_half_ratio(numerator: float, denominator: float) -> float:
    return (numerator + 0.5) / (denominator + 0.5)


@dataclass
class ArchitectureObject:
    species: str
    object_type: str
    object_id: str
    segments: List[Segment]


print('=' * 72)
print('PAIRWISE ARCHITECTURE TEST v2')
print('=' * 72)

mat_path = first_existing(MAT_CANDIDATES, 'MAT')
master3_path = first_existing(MASTER3_CANDIDATES, 'alignment3 master')
mat = pd.read_csv(mat_path)
master3 = pd.read_csv(master3_path)

required_cols = {'species', 'object_type', 'alignment3_start', 'alignment3_end'}
missing = required_cols - set(mat.columns)
if missing:
    raise ValueError(f'В MAT отсутствуют столбцы: {sorted(missing)}')
if 'alignment3_column' not in master3.columns:
    raise ValueError('В alignment3 master отсутствует столбец alignment3_column')

alignment3_length = int(pd.to_numeric(master3['alignment3_column'], errors='raise').max())
print('MAT:', mat_path)
print('alignment3 master:', master3_path)
print('Длина alignment3:', alignment3_length)


def extract_interval_objects(data: pd.DataFrame, species: str, object_type: str, prefix: str) -> List[ArchitectureObject]:
    layer = data[(data['species'] == species) & (data['object_type'] == object_type)].copy()
    starts = pd.to_numeric(layer['alignment3_start'], errors='coerce')
    ends = pd.to_numeric(layer['alignment3_end'], errors='coerce')
    valid = starts.notna() & ends.notna()
    objects: List[ArchitectureObject] = []
    for counter, (start, end) in enumerate(zip(starts[valid], ends[valid]), start=1):
        objects.append(ArchitectureObject(species, object_type, f'{prefix}_{species}_{counter:03d}', [normalize_segment(int(start), int(end))]))
    return objects


def extract_cleavage_objects(data: pd.DataFrame, species: str) -> List[ArchitectureObject]:
    layer = data[(data['species'] == species) & (data['object_type'] == 'CLEAVAGE_SITE')].copy()
    left = pd.to_numeric(layer['alignment3_boundary_left'], errors='coerce') if 'alignment3_boundary_left' in layer.columns else pd.Series(np.nan, index=layer.index)
    right = pd.to_numeric(layer['alignment3_boundary_right'], errors='coerce') if 'alignment3_boundary_right' in layer.columns else pd.Series(np.nan, index=layer.index)
    left = left.fillna(pd.to_numeric(layer['alignment3_start'], errors='coerce'))
    right = right.fillna(pd.to_numeric(layer['alignment3_end'], errors='coerce'))
    valid = left.notna() & right.notna()
    objects: List[ArchitectureObject] = []
    for counter, (start, end) in enumerate(zip(left[valid], right[valid]), start=1):
        objects.append(ArchitectureObject(species, 'CLEAVAGE_SITE', f'CLEAVAGE_{species}_{counter:03d}', [normalize_segment(int(start), int(end))]))
    return objects


def extract_structure_objects(data: pd.DataFrame, species: str, model: str) -> List[ArchitectureObject]:
    layer = data[
        (data['species'] == species)
        & (data['object_type'] == 'RNA_STRUCTURE')
        & data['source_column'].astype(str).str.lower().str.contains(model.lower(), na=False)
    ].copy()
    starts = pd.to_numeric(layer['alignment3_start'], errors='coerce')
    ends = pd.to_numeric(layer['alignment3_end'], errors='coerce')
    valid = starts.notna() & ends.notna()
    objects: List[ArchitectureObject] = []
    for counter, (start, end) in enumerate(zip(starts[valid], ends[valid]), start=1):
        objects.append(ArchitectureObject(species, f'RNA_STRUCTURE_{model.upper()}', f'{model.upper()}_{species}_{counter:04d}', [normalize_segment(int(start), int(end))]))
    return objects


def objects_to_intervals(objects: Sequence[ArchitectureObject]) -> np.ndarray:
    intervals: List[Segment] = []
    for obj in objects:
        intervals.extend(obj.segments)
    return np.asarray(intervals, dtype=int) if intervals else np.empty((0, 2), dtype=int)


def objects_to_boundaries(objects: Sequence[ArchitectureObject]) -> np.ndarray:
    boundaries = set()
    for obj in objects:
        for start, end in obj.segments:
            boundaries.add(int(start)); boundaries.add(int(end))
    return np.asarray(sorted(boundaries), dtype=int)


ecm_objects: Dict[str, List[ArchitectureObject]] = {}
apa_objects: Dict[str, List[ArchitectureObject]] = {}
mfe_objects: Dict[str, List[ArchitectureObject]] = {}
centroid_objects: Dict[str, List[ArchitectureObject]] = {}

for species in SPECIES_ORDER:
    ecm_objects[species] = extract_interval_objects(mat, species, 'ECM_BLOCK', 'ECM')
    current_apa: List[ArchitectureObject] = []
    current_apa += extract_interval_objects(mat, species, 'UGUA', 'UGUA')
    current_apa += extract_interval_objects(mat, species, 'PAS', 'PAS')
    current_apa += extract_cleavage_objects(mat, species)
    current_apa += extract_interval_objects(mat, species, 'DOWNSTREAM_ELEMENT', 'DSE')
    apa_objects[species] = current_apa
    mfe_objects[species] = extract_structure_objects(mat, species, 'mfe')
    centroid_objects[species] = extract_structure_objects(mat, species, 'centroid')

structure_intervals = {'ANY': {}, 'MFE': {}, 'CENTROID': {}}
structure_boundaries = {'ANY': {}, 'MFE': {}, 'CENTROID': {}}
for species in SPECIES_ORDER:
    structure_intervals['MFE'][species] = objects_to_intervals(mfe_objects[species])
    structure_intervals['CENTROID'][species] = objects_to_intervals(centroid_objects[species])
    structure_intervals['ANY'][species] = objects_to_intervals(mfe_objects[species] + centroid_objects[species])
    structure_boundaries['MFE'][species] = objects_to_boundaries(mfe_objects[species])
    structure_boundaries['CENTROID'][species] = objects_to_boundaries(centroid_objects[species])
    structure_boundaries['ANY'][species] = objects_to_boundaries(mfe_objects[species] + centroid_objects[species])

inventory_rows = []
for species in SPECIES_ORDER:
    inventory_rows.append({
        'species': species,
        'n_ECM': len(ecm_objects[species]),
        'n_APA': len(apa_objects[species]),
        'n_UGUA': sum(o.object_type == 'UGUA' for o in apa_objects[species]),
        'n_PAS': sum(o.object_type == 'PAS' for o in apa_objects[species]),
        'n_CLEAVAGE': sum(o.object_type == 'CLEAVAGE_SITE' for o in apa_objects[species]),
        'n_DSE': sum(o.object_type == 'DOWNSTREAM_ELEMENT' for o in apa_objects[species]),
        'n_MFE_structure': len(mfe_objects[species]),
        'n_centroid_structure': len(centroid_objects[species]),
    })
inventory = pd.DataFrame(inventory_rows)
print('\nИнвентаризация объектов:')
print(inventory.to_string(index=False))


def relation_to_structure(segments: Sequence[Segment], species: str, model: str) -> Dict[str, int]:
    intervals = structure_intervals[model][species]
    boundaries = structure_boundaries[model][species]
    overlap = interval_overlap(segments, intervals)
    boundary_distance = minimum_boundary_distance(segments, boundaries)
    return {
        'overlap': int(overlap),
        'boundary': int(boundary_distance <= BOUNDARY_WINDOW),
        'center': int(overlap and boundary_distance > BOUNDARY_WINDOW),
    }

metric_meta: List[Dict[str, str]] = []

def register(metric: str, edge_id: str, family: str, interpretation: str):
    metric_meta.append({'metric': metric, 'edge_id': edge_id, 'metric_family': family, 'interpretation': interpretation})

register('ECM_STRUCTURE__ANY__boundary_to_center', 'ECM__RNA_STRUCTURE', 'primary', 'ECM boundary/center ratio relative to any predicted structure.')
register('ECM_STRUCTURE__ANY__overlap', 'ECM__RNA_STRUCTURE', 'primary', 'ECM overlapping any predicted structural element.')
register('ECM_APA__all__overlap', 'ECM__APA', 'primary', 'APA objects overlapping ECM.')
register('ECM_APA__all__within_10', 'ECM__APA', 'primary', 'APA objects within 10 alignment columns of ECM.')
register('STRUCTURE_APA__ANY__boundary_to_center', 'RNA_STRUCTURE__APA', 'primary', 'APA boundary/center ratio relative to any predicted structure.')
register('STRUCTURE_APA__ANY__overlap', 'RNA_STRUCTURE__APA', 'primary', 'APA objects overlapping any predicted structural element.')

for model in ['ANY', 'MFE', 'CENTROID']:
    for rel in ['boundary', 'center']:
        register(f'ECM_STRUCTURE__{model}__{rel}', 'ECM__RNA_STRUCTURE', 'exploratory', f'ECM {rel} relative to {model} structure.')
    if model != 'ANY':
        register(f'ECM_STRUCTURE__{model}__overlap', 'ECM__RNA_STRUCTURE', 'exploratory', f'ECM overlap with {model} structure.')
        register(f'ECM_STRUCTURE__{model}__boundary_to_center', 'ECM__RNA_STRUCTURE', 'exploratory', f'ECM boundary/center ratio for {model}.')
register('ECM_APA__all__within_5', 'ECM__APA', 'exploratory', 'APA objects within 5 columns of ECM.')
for apa_type in APA_TYPES:
    for rel in ['overlap', 'within_10']:
        register(f'ECM_APA__{apa_type}__{rel}', 'ECM__APA', 'exploratory', f'{apa_type} {rel} relative to ECM.')
for model in ['ANY', 'MFE', 'CENTROID']:
    for rel in ['boundary', 'center']:
        register(f'STRUCTURE_APA__{model}__{rel}', 'RNA_STRUCTURE__APA', 'exploratory', f'APA {rel} relative to {model} structure.')
    if model != 'ANY':
        register(f'STRUCTURE_APA__{model}__overlap', 'RNA_STRUCTURE__APA', 'exploratory', f'APA overlap with {model} structure.')
        register(f'STRUCTURE_APA__{model}__boundary_to_center', 'RNA_STRUCTURE__APA', 'exploratory', f'APA boundary/center ratio for {model}.')
for apa_type in APA_TYPES:
    for rel in ['boundary', 'center']:
        register(f'STRUCTURE_APA__{apa_type}__{rel}', 'RNA_STRUCTURE__APA', 'exploratory', f'{apa_type} {rel} relative to any structure.')

metric_metadata = pd.DataFrame(metric_meta)
pairwise_metrics = metric_metadata['metric'].tolist()
if metric_metadata['metric'].duplicated().any():
    raise ValueError('Есть дублированные метрики')


def count_pairwise_metrics(moving_ecm: Dict[str, List[ArchitectureObject]], moving_apa: Dict[str, List[ArchitectureObject]]) -> Dict[str, float]:
    result = {metric: 0.0 for metric in pairwise_metrics}

    ecm_structure = {m: {'overlap': 0, 'boundary': 0, 'center': 0} for m in ['ANY', 'MFE', 'CENTROID']}
    for species in SPECIES_ORDER:
        for ecm in moving_ecm[species]:
            for model in ['ANY', 'MFE', 'CENTROID']:
                rel = relation_to_structure(ecm.segments, species, model)
                for key in ['overlap', 'boundary', 'center']:
                    ecm_structure[model][key] += rel[key]
    result['ECM_STRUCTURE__ANY__overlap'] = ecm_structure['ANY']['overlap']
    result['ECM_STRUCTURE__ANY__boundary_to_center'] = add_half_ratio(ecm_structure['ANY']['boundary'], ecm_structure['ANY']['center'])
    for model in ['ANY', 'MFE', 'CENTROID']:
        result[f'ECM_STRUCTURE__{model}__boundary'] = ecm_structure[model]['boundary']
        result[f'ECM_STRUCTURE__{model}__center'] = ecm_structure[model]['center']
        if model != 'ANY':
            result[f'ECM_STRUCTURE__{model}__overlap'] = ecm_structure[model]['overlap']
            result[f'ECM_STRUCTURE__{model}__boundary_to_center'] = add_half_ratio(ecm_structure[model]['boundary'], ecm_structure[model]['center'])

    ecm_apa = {'all': {'overlap': 0, 'within_5': 0, 'within_10': 0}}
    for t in APA_TYPES:
        ecm_apa[t] = {'overlap': 0, 'within_10': 0}
    for species in SPECIES_ORDER:
        ecm_intervals = objects_to_intervals(ecm_objects[species])
        for apa in moving_apa[species]:
            overlap = interval_overlap(apa.segments, ecm_intervals)
            dist = minimum_interval_distance(apa.segments, ecm_intervals)
            ecm_apa['all']['overlap'] += int(overlap)
            ecm_apa['all']['within_5'] += int(dist <= 5)
            ecm_apa['all']['within_10'] += int(dist <= 10)
            ecm_apa[apa.object_type]['overlap'] += int(overlap)
            ecm_apa[apa.object_type]['within_10'] += int(dist <= 10)
    result['ECM_APA__all__overlap'] = ecm_apa['all']['overlap']
    result['ECM_APA__all__within_5'] = ecm_apa['all']['within_5']
    result['ECM_APA__all__within_10'] = ecm_apa['all']['within_10']
    for t in APA_TYPES:
        result[f'ECM_APA__{t}__overlap'] = ecm_apa[t]['overlap']
        result[f'ECM_APA__{t}__within_10'] = ecm_apa[t]['within_10']

    structure_apa = {m: {'overlap': 0, 'boundary': 0, 'center': 0} for m in ['ANY', 'MFE', 'CENTROID']}
    structure_apa_type = {t: {'boundary': 0, 'center': 0} for t in APA_TYPES}
    for species in SPECIES_ORDER:
        for apa in moving_apa[species]:
            for model in ['ANY', 'MFE', 'CENTROID']:
                rel = relation_to_structure(apa.segments, species, model)
                for key in ['overlap', 'boundary', 'center']:
                    structure_apa[model][key] += rel[key]
                if model == 'ANY':
                    structure_apa_type[apa.object_type]['boundary'] += rel['boundary']
                    structure_apa_type[apa.object_type]['center'] += rel['center']
    result['STRUCTURE_APA__ANY__overlap'] = structure_apa['ANY']['overlap']
    result['STRUCTURE_APA__ANY__boundary_to_center'] = add_half_ratio(structure_apa['ANY']['boundary'], structure_apa['ANY']['center'])
    for model in ['ANY', 'MFE', 'CENTROID']:
        result[f'STRUCTURE_APA__{model}__boundary'] = structure_apa[model]['boundary']
        result[f'STRUCTURE_APA__{model}__center'] = structure_apa[model]['center']
        if model != 'ANY':
            result[f'STRUCTURE_APA__{model}__overlap'] = structure_apa[model]['overlap']
            result[f'STRUCTURE_APA__{model}__boundary_to_center'] = add_half_ratio(structure_apa[model]['boundary'], structure_apa[model]['center'])
    for t in APA_TYPES:
        result[f'STRUCTURE_APA__{t}__boundary'] = structure_apa_type[t]['boundary']
        result[f'STRUCTURE_APA__{t}__center'] = structure_apa_type[t]['center']
    return result

observed_pairwise = count_pairwise_metrics(ecm_objects, apa_objects)
observed_table = metric_metadata.copy()
observed_table['observed'] = [observed_pairwise[m] for m in pairwise_metrics]
print('\nНаблюдаемые метрики:')
print(observed_table.to_string(index=False))


def shift_collection(objects_by_species: Dict[str, List[ArchitectureObject]], shift: int) -> Dict[str, List[ArchitectureObject]]:
    out: Dict[str, List[ArchitectureObject]] = {}
    for species in SPECIES_ORDER:
        out[species] = [
            ArchitectureObject(o.species, o.object_type, o.object_id, shift_segments_circular(o.segments, shift, alignment3_length))
            for o in objects_by_species[species]
        ]
    return out

pairwise_shifts = np.arange(1, alignment3_length, dtype=int)
metric_to_index = {m: i for i, m in enumerate(pairwise_metrics)}
null_pairwise = np.zeros((len(pairwise_shifts), len(pairwise_metrics)), dtype=float)
metric_groups = {
    'ECM_STRUCTURE': [m for m in pairwise_metrics if m.startswith('ECM_STRUCTURE__')],
    'ECM_APA': [m for m in pairwise_metrics if m.startswith('ECM_APA__')],
    'STRUCTURE_APA': [m for m in pairwise_metrics if m.startswith('STRUCTURE_APA__')],
}

print('\nНачинается полный перебор циклических сдвигов...')
start_time = time.time()
for shift_index, shift in enumerate(pairwise_shifts):
    shifted_ecm = shift_collection(ecm_objects, int(shift))
    shifted_apa = shift_collection(apa_objects, int(shift))
    ecm_structure_values = count_pairwise_metrics(shifted_ecm, apa_objects)
    apa_shift_values = count_pairwise_metrics(ecm_objects, shifted_apa)
    for m in metric_groups['ECM_STRUCTURE']:
        null_pairwise[shift_index, metric_to_index[m]] = ecm_structure_values[m]
    for m in metric_groups['ECM_APA'] + metric_groups['STRUCTURE_APA']:
        null_pairwise[shift_index, metric_to_index[m]] = apa_shift_values[m]
    if (shift_index + 1) % 100 == 0 or (shift_index + 1) == len(pairwise_shifts):
        print(f'Выполнено: {shift_index + 1}/{len(pairwise_shifts)}; {time.time() - start_time:.1f} сек.')

result_rows = []
for i, metric in enumerate(pairwise_metrics):
    meta = metric_metadata.loc[metric_metadata['metric'] == metric].iloc[0]
    observed = float(observed_pairwise[metric])
    null_values = null_pairwise[:, i]
    null_mean = float(np.mean(null_values))
    null_sd = float(np.std(null_values, ddof=1))
    p_enrichment = (1 + int(np.sum(null_values >= observed))) / (len(null_values) + 1)
    p_depletion = (1 + int(np.sum(null_values <= observed))) / (len(null_values) + 1)
    fold = observed / null_mean if null_mean > 0 else np.inf
    z = (observed - null_mean) / null_sd if null_sd > 0 else np.nan
    result_rows.append({
        'edge_id': meta['edge_id'], 'metric': metric, 'metric_family': meta['metric_family'], 'interpretation': meta['interpretation'],
        'observed': observed, 'null_mean': null_mean, 'null_median': float(np.median(null_values)), 'null_sd': null_sd,
        'null_q025': float(np.quantile(null_values, 0.025)), 'null_q975': float(np.quantile(null_values, 0.975)),
        'fold_observed_expected': fold, 'z_score': z, 'p_enrichment': p_enrichment, 'p_depletion': p_depletion,
        'n_shifts': len(null_values),
    })

pairwise_results = pd.DataFrame(result_rows)
pairwise_results['q_enrichment_BH'] = np.nan
pairwise_results['q_depletion_BH'] = np.nan
for (edge_id, family), idx in pairwise_results.groupby(['edge_id', 'metric_family']).groups.items():
    idx = list(idx)
    pairwise_results.loc[idx, 'q_enrichment_BH'] = benjamini_hochberg(pairwise_results.loc[idx, 'p_enrichment'].to_numpy())
    pairwise_results.loc[idx, 'q_depletion_BH'] = benjamini_hochberg(pairwise_results.loc[idx, 'p_depletion'].to_numpy())


def classify(row):
    if row['q_enrichment_BH'] < 0.05 and row['fold_observed_expected'] > 1:
        return 'significant_enrichment'
    if row['q_depletion_BH'] < 0.05 and row['fold_observed_expected'] < 1:
        return 'significant_depletion'
    if row['p_enrichment'] < 0.05 and row['fold_observed_expected'] > 1:
        return 'nominal_enrichment'
    if row['p_depletion'] < 0.05 and row['fold_observed_expected'] < 1:
        return 'nominal_depletion'
    return 'not_significant'

pairwise_results['statistical_class'] = pairwise_results.apply(classify, axis=1)

edge_summary_rows = []
for edge_id, subset in pairwise_results.groupby('edge_id'):
    primary = subset[subset['metric_family'] == 'primary'].copy()
    sig_primary = primary[primary['statistical_class'].isin(['significant_enrichment', 'significant_depletion'])]
    nom_primary = primary[primary['statistical_class'].isin(['nominal_enrichment', 'nominal_depletion'])]
    sig_any = subset[subset['statistical_class'].isin(['significant_enrichment', 'significant_depletion'])]
    nom_any = subset[subset['statistical_class'].isin(['nominal_enrichment', 'nominal_depletion'])]
    if not sig_primary.empty:
        grade, status = 'A', 'FDR_significant_primary_relation'
    elif not nom_primary.empty:
        grade, status = 'B', 'nominal_primary_relation'
    elif not sig_any.empty:
        grade, status = 'C', 'FDR_significant_exploratory_relation'
    elif not nom_any.empty:
        grade, status = 'C', 'nominal_exploratory_relation'
    else:
        grade, status = 'D', 'tested_no_pairwise_association'
    candidates = primary if not primary.empty else subset.copy()
    candidates = candidates.copy()
    candidates['best_q'] = candidates[['q_enrichment_BH', 'q_depletion_BH']].min(axis=1)
    candidates['best_p'] = candidates[['p_enrichment', 'p_depletion']].min(axis=1)
    strongest = candidates.sort_values(['best_q', 'best_p', 'metric']).iloc[0]
    edge_summary_rows.append({
        'edge_id': edge_id, 'evidence_grade': grade, 'evidence_status': status,
        'strongest_metric': strongest['metric'], 'strongest_metric_family': strongest['metric_family'],
        'observed': strongest['observed'], 'null_mean': strongest['null_mean'],
        'effect_ratio': strongest['fold_observed_expected'], 'p_enrichment': strongest['p_enrichment'],
        'p_depletion': strongest['p_depletion'], 'q_enrichment': strongest['q_enrichment_BH'],
        'q_depletion': strongest['q_depletion_BH'], 'statistical_class': strongest['statistical_class'],
    })
pairwise_edge_summary = pd.DataFrame(edge_summary_rows)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
BACKUP_DIR.mkdir(parents=True, exist_ok=True)
BAN_V2_BACKUP_DIR.mkdir(parents=True, exist_ok=True)
results_csv = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2_results.csv'
edge_summary_csv = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2_edge_summary.csv'
inventory_csv = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2_inventory.csv'
results_xlsx = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2.xlsx'
null_npz = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2_null.npz'
qc_txt = OUTPUT_DIR / 'Ras85D_Pairwise_Architecture_Test_v2_QC.txt'

pairwise_results.to_csv(results_csv, index=False)
pairwise_edge_summary.to_csv(edge_summary_csv, index=False)
inventory.to_csv(inventory_csv, index=False)
np.savez_compressed(null_npz, metrics=np.asarray(pairwise_metrics, dtype=object), shifts=pairwise_shifts, null=null_pairwise, alignment3_length=np.asarray([alignment3_length]))

with pd.ExcelWriter(results_xlsx, engine='openpyxl') as writer:
    pairwise_results.to_excel(writer, sheet_name='all_tests', index=False)
    pairwise_edge_summary.to_excel(writer, sheet_name='edge_summary', index=False)
    inventory.to_excel(writer, sheet_name='object_inventory', index=False)
    observed_table.to_excel(writer, sheet_name='observed_metrics', index=False)
    for edge_id, subset in pairwise_results.groupby('edge_id'):
        subset.to_excel(writer, sheet_name=edge_id.replace('__', '_')[:31], index=False)

qc_txt.write_text('\n'.join([
    'Ras85D Pairwise Architecture Test v2', '=' * 72,
    f'MAT: {mat_path}', f'alignment3 master: {master3_path}',
    f'alignment3 length: {alignment3_length}', f'non-zero circular shifts: {len(pairwise_shifts)}',
    f'tested metrics: {len(pairwise_metrics)}', '',
    'FDR was applied separately within each BAN edge and metric family.',
    f'Boundary association: distance <= {BOUNDARY_WINDOW} alignment columns.', '',
    'Edge summary:', pairwise_edge_summary.to_string(index=False), '',
    'Limitation: the test evaluates relations relative to the observed architectural landscape and does not imply direct 3D contact.'
]), encoding='utf-8')

ban_v2_files: List[Path] = []
ban_edge_path = next((p for p in BAN_EDGE_CANDIDATES if p.exists()), None)
ban_node_path = next((p for p in BAN_NODE_CANDIDATES if p.exists()), None)
if ban_edge_path is not None:
    ban_edges_v2 = pd.read_csv(ban_edge_path)
    for _, s in pairwise_edge_summary.iterrows():
        mask = ban_edges_v2['edge_id'].eq(s['edge_id'])
        if mask.any():
            ban_edges_v2.loc[mask, 'evidence_grade'] = s['evidence_grade']
            ban_edges_v2.loc[mask, 'evidence_status'] = s['evidence_status']
            ban_edges_v2.loc[mask, 'effect_direction'] = 'positive' if s['effect_ratio'] > 1 else 'negative' if s['effect_ratio'] < 1 else 'none_detected'
            ban_edges_v2.loc[mask, 'effect_summary'] = f"{s['strongest_metric']}: O/E={s['effect_ratio']:.3f}"
            ban_edges_v2.loc[mask, 'statistical_summary'] = f"{s['statistical_class']}; p_enr={s['p_enrichment']:.4g}; p_dep={s['p_depletion']:.4g}; q_enr={s['q_enrichment']:.4g}; q_dep={s['q_depletion']:.4g}."
            ban_edges_v2.loc[mask, 'primary_result'] = 'Pairwise Architecture Test v2 quantified this relation using exhaustive configuration-preserving circular shifts.'
    ban_edges_v2_csv = OUTPUT_DIR / 'Ras85D_BAN_v2_edges.csv'
    ban_edges_v2.to_csv(ban_edges_v2_csv, index=False)
    ban_v2_files.append(ban_edges_v2_csv)
else:
    ban_edges_v2 = pd.DataFrame()

if ban_node_path is not None:
    ban_nodes_v2 = pd.read_csv(ban_node_path)
    ban_nodes_v2_csv = OUTPUT_DIR / 'Ras85D_BAN_v2_nodes.csv'
    ban_nodes_v2.to_csv(ban_nodes_v2_csv, index=False)
    ban_v2_files.append(ban_nodes_v2_csv)
else:
    ban_nodes_v2 = pd.DataFrame()

if not ban_edges_v2.empty:
    grade_score = {'A': 4, 'B': 3, 'C': 2, 'D': 1, 'U': 0}
    node_order = ['IEL', 'ECM', 'RNA_STRUCTURE', 'APA']
    adjacency = pd.DataFrame(0, index=node_order, columns=node_order, dtype=int)
    for _, edge in ban_edges_v2.iterrows():
        s, t = edge['source_node'], edge['target_node']
        if s in adjacency.index and t in adjacency.columns:
            score = grade_score.get(edge['evidence_grade'], 0)
            adjacency.loc[s, t] = score
            adjacency.loc[t, s] = score
    adjacency_csv = OUTPUT_DIR / 'Ras85D_BAN_v2_adjacency.csv'
    adjacency.to_csv(adjacency_csv)
    ban_v2_files.append(adjacency_csv)

main_files = [results_csv, edge_summary_csv, inventory_csv, results_xlsx, null_npz, qc_txt]
for source in main_files:
    target = BACKUP_DIR / source.name
    shutil.copy2(source, target)
    print('Скопирован:', target)
for source in ban_v2_files:
    target = BAN_V2_BACKUP_DIR / source.name
    shutil.copy2(source, target)
    print('Скопирован:', target)

print('\n' + '=' * 72)
print('PAIRWISE ARCHITECTURE TEST v2 ЗАВЕРШЁН')
print('=' * 72)
print('\nИтоги по рёбрам BAN:')
print(pairwise_edge_summary.to_string(index=False))
print('\nЗначимые и номинальные метрики:')
interesting = pairwise_results[pairwise_results['statistical_class'] != 'not_significant']
print('Нет значимых или номинальных результатов.' if interesting.empty else interesting.to_string(index=False))
print('\nBackup:')
print(BACKUP_DIR)
print(BAN_V2_BACKUP_DIR)
