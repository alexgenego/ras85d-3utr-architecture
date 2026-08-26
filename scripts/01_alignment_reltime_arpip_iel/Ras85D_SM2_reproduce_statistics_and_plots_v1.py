# Reproduce Table S2 statistics and Figure S1 panels B-E.
# Input: Ras85D_Table_S2_tip_rates_v1.csv
# Requires: Python 3.11+, numpy, scipy, matplotlib
import csv
from pathlib import Path
import numpy as np
from scipy import stats

path = Path("Ras85D_Table_S2_tip_rates_v1.csv")
with path.open(encoding="utf-8") as f:
    rows = list(csv.DictReader(f))

partial_rel = np.array([float(r["partial_MEGA_RelTime_relative_rate"]) for r in rows])
complete_rel = np.array([float(r["complete_MEGA_RelTime_relative_rate"]) for r in rows])
partial_cum = np.array([float(r["partial_cumulative_root_to_tip_rate_subs_per_site_per_year"]) for r in rows])
complete_cum = np.array([float(r["complete_cumulative_root_to_tip_rate_subs_per_site_per_year"]) for r in rows])

print("Terminal relative rate")
print("Wilcoxon:", stats.wilcoxon(partial_rel, complete_rel))
print("Pearson:", stats.pearsonr(partial_rel, complete_rel))
print("Spearman:", stats.spearmanr(partial_rel, complete_rel))

print("\nCumulative root-to-tip rate")
print("Wilcoxon:", stats.wilcoxon(partial_cum, complete_cum))
print("Pearson:", stats.pearsonr(partial_cum, complete_cum))
print("Spearman:", stats.spearmanr(partial_cum, complete_cum))
