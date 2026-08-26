# PAT, AEM, Coverage Enrichment and Functional Distance

This directory contains the manuscript-associated implementations and reproducibility wrappers for the Pairwise Architecture Test (PAT), Architecture Evidence Matrix (AEM), Coverage Enrichment Test, and Functional Distance Test.

## Recommended PAT/AEM reproduction route

Use:

- `Ras85D_prepare_authoritative_PAT_AEM_input_v2.py`;
- `Ras85D_run_authoritative_PAT_AEM_v2.py`;
- the original statistical implementations `Pairwise_Architecture_Test_v2.py` and `build_ADT_v1.py`.

The v2 preparation/wrapper pair corrects packaging/compatibility issues in the archived Supplementary wrappers only. It adds the compatibility fields required by the original PAT/AEM scripts, redirects historical Google Colab filesystem paths to the supplied workspace, and validates the final total of 43 PAT metrics. The statistical calculations in `Pairwise_Architecture_Test_v2.py` and `build_ADT_v1.py` are not altered.

With the manuscript-associated final registry and the documented reconstructed IEL inputs, this route reproduces the preserved PAT results (43 metrics and three edge summaries) and the AEM output (102 IEL × species rows and 34 consensus rows). The rerun was checked against the preserved final machine-readable outputs.

The files `Ras85D_prepare_authoritative_PAT_AEM_input_v1.py` and `Ras85D_run_authoritative_PAT_AEM_v1.py` are retained as the original archived Supplementary wrappers for provenance. They should not be used as the final public rerun route because the v1 preparation omitted compatibility columns and the v1 wrapper expected 41 rather than the final 43 PAT metrics.

`Ras85D_SM9_validate_summary_and_reconciliation_v1.py` and `Ras85D_SM8_preflight_and_validate_v1.py` are historical Supplementary-assembly/QC helpers. They document earlier reconciliation stages and are not the validators for the final 43-metric PAT/AEM result set.

## Coverage Enrichment and Functional Distance

`Ras85D_run_Table_S9_authoritative_v1.py` is the recommended manuscript-associated runner for the final Coverage Enrichment and Functional Distance analyses. It uses the final 510-object *D. melanogaster* registry, 1-based inclusive 3′ UTR coordinates 1–964, 10,000 interval-preserving permutations, and seed 20260718.

A physical rerun reproduced the manuscript-associated numerical statistics. The six FDR-significant Coverage rows and all 12 Functional Distance rows match the preserved final output tables.

Historical filenames containing `authoritative` are retained for provenance; the term is not used as scientific terminology in the manuscript.
