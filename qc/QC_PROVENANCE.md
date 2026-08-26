# Quality-control provenance

This directory contains manuscript-relevant validation records used to audit the public Ras85D 3′ UTR reproducibility package.

## Scope

The archived QC records cover:

- reconstruction of the RelTime tip-rate CSV from the preserved Supplementary Table S2 worksheet and reproduction of the reported paired statistics;
- reproduction of the 34-IEL spatial summary;
- validation of ECM, mobile-element, microsatellite, and functional-object registries;
- validation of preserved Boundary Pair Universe (BPU) summary counts;
- complete PAT/AEM rerun from the recovered minimal inputs, followed by direct comparison with preserved final outputs;
- independent rerun of Coverage Enrichment Test v2 and Functional Distance Test v1;
- exact reproduction of the manuscript-associated five-state multiloop dispersion randomization after preserving the original first-appearance species order.

## PAT/AEM recovery audit

Two historical intermediate PAT/AEM input files were reconstructed from preserved final IEL/projection evidence. The authoritative rerun produced 43 PAT metrics, 3 PAT edges, 102 AEM long-form rows, and 34 AEM consensus rows. Direct comparison with the preserved final PAT and AEM tables showed zero numeric and text discrepancies in the compared final outputs.

## Boundary Pair Universe limitation

The historical normalized 525-object BPU architecture input and observed-interface catalogue were not retained after working-archive cleanup. They are not reconstructed or substituted here. The public archive preserves the final BPU scripts, historical SHA-256 checksums recorded for the two missing inputs, the manuscript-associated BPU QC summary, and final compact boundary-neighborhood outputs.

## Multiloop dispersion script correction

The archived Supplementary Script S5 v1 used `np.unique()` to enumerate species, which changed the species iteration order relative to the manuscript-associated run and therefore changed the deterministic random-number stream despite using the same seed. Public reproducibility Script S5 v2 preserves species in first-appearance order. With that compatibility correction, the primary five-state statistics and primary contrast reproduce the preserved manuscript values exactly (within floating-point representation).

## Interpretation

These QC records document computational reproducibility and consistency of the archived analysis products. They do not add new biological analyses beyond those reported in the manuscript.
