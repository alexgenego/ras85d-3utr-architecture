# Multiloop-context dispersion

Final script for the primary five-state between-context dispersion analysis.

Typical execution:

```bash
python Ras85D_run_multiloop_dispersion_v1.py \
  --usage 29A_3_08B_Tissue_Resolved_APA_Cleavage_Efficiency_v2.xlsx \
  --structure 29A_3_09B_multiloop_state_data.xlsx \
  --outdir dispersion_output \
  --permutations 100000 \
  --seed 20260727
```

Python dependencies include NumPy, pandas, and openpyxl.
