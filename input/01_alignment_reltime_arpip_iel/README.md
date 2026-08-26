# Alignment, RelTime, ARPIP and IEL inputs

This directory contains preserved manuscript-associated sequence/alignment inputs, corrected tree/name records, the compact IEL datasets, and the terminal RelTime vectors used by the public reproducibility scripts.

The file `Ras85D_Table_S2_tip_rates_v1.csv` was reconstructed losslessly from the preserved `Tip_rates` worksheet of Supplementary Table S2 because the standalone CSV was not present in the compact Supplementary archive. The four rate vectors used by the script were validated by rerunning `Ras85D_SM2_reproduce_statistics_and_plots_v1.py`; the reported Wilcoxon, Pearson and Spearman statistics reproduced the final Supplementary values.

The ARPIP program itself and large exploratory/intermediate ARPIP work files are not redistributed here. Final branch/IEL outputs are archived separately under `output/`.
