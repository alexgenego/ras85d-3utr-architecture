# Multiloop-context dispersion

This directory contains the compact manuscript-associated implementation of the primary five-state between-context dispersion analysis.

## Recommended primary reproduction route

Use `Ras85D_run_multiloop_dispersion_v2.py`.

The v2 script differs from the archived Supplementary `Ras85D_run_multiloop_dispersion_v1.py` only in how species blocks are ordered when constructing the fixed randomization stream. The manuscript-associated analysis used the first-appearance source-table order:

1. `D.melanogaster`;
2. `D.yakuba`;
3. `D.virilis`.

The archived v1 script used `np.unique(species)`, which sorts species alphabetically and therefore changes the sequence of random draws despite the same seed. No statistical model, statistic, weighting rule, permutation constraint, or input data were changed in v2.

Typical execution:

```bash
python Ras85D_run_multiloop_dispersion_v2.py \
  --usage 29A_3_08B_Tissue_Resolved_APA_Cleavage_Efficiency_v2.xlsx \
  --structure 29A_3_09B_multiloop_state_data.xlsx \
  --outdir dispersion_output \
  --permutations 100000 \
  --seed 20260727
```

With NumPy 2.3.5 and the archived inputs, v2 reproduces the preserved primary five-state result table and the prespecified `DIFFERENT_MULTILOOP_ARMS − SAME_MULTILOOP_INTERVAL` all-context contrast exactly, including permutation p-values and FDR values.

The archived v1 script is retained unchanged for Supplementary provenance but is not the recommended exact-rerun implementation.

The compact v2 script reproduces the primary five-state/all-context analysis. The manuscript's additional tissue/embryo, model-specific, pair-scope, species-specific, and leave-one-species-out sensitivity outputs are preserved in `output/`; this compact script does not regenerate every sensitivity analysis.

Python dependencies include NumPy, pandas, and openpyxl.
