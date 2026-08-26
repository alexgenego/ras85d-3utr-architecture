# Input provenance and recovery status

This public input set contains the final manuscript-relevant files needed by the archived custom scripts, while avoiding duplication of large intermediate working archives.

## Preserved original inputs

The sequence/alignment files, final harmonized functional-object registry, compact ECM/mobile-element/microsatellite catalogues, IEL tables, coordinate rules, and the two multiloop-dispersion workbooks are copied from the final Supplementary Data packages assembled for the manuscript.

## RelTime vector recovery

The standalone `Ras85D_Table_S2_tip_rates_v1.csv` was absent from the compact Supplementary Data archive, although the same 36-taxon vectors were retained in Supplementary Table S2. The CSV in this repository was reconstructed directly from the preserved `Tip_rates` worksheet and converted to the column names expected by the archived reproducibility script. Rerunning the script reproduced the final paired Wilcoxon, Pearson and Spearman statistics.

## PAT/AEM recovery

Two historical intermediate files required by the original PAT/AEM wrapper were lost during archive cleanup. Minimal recovery files were reconstructed from preserved final IEL/projection evidence and then tested by a complete PAT/AEM rerun. The recovered inputs reproduce the preserved final PAT and AEM results, and are explicitly marked `reconstructed` in their filenames.

## Boundary Pair Universe limitation

The historical normalized 525-object BPU architecture input and observed-interface catalogue are no longer available. They have not been approximated using a later registry. Their recorded SHA-256 values, final scripts, and preserved QC/result evidence are retained so that the limitation is transparent.
