# PAT/AEM recovery inputs

The two historical intermediate PAT/AEM input CSV files were not retained after project-archive cleanup. They were reconstructed from preserved final manuscript-associated tables without changing the biological annotations used by the final analysis:

- `Ras85D_IEL_master_table_reconstructed_minimal.csv`
- `arpip_34_regions_alignment3_species_specific_reconstructed.csv`

The final object registry used by PAT/AEM and by Coverage Enrichment / Functional Distance is stored once at:

`../02_feature_catalogue_registry_qc/29A_3_01_final_object_registry.csv`

Recovery was validated by rerunning the original PAT and AEM statistical logic. The rerun produced 43 PAT metrics and 3 tested edges, with numerical PAT results matching the preserved final table to better than 1e-12, and a 34-row AEM consensus table matching the preserved final consensus across all 37 columns.

These reconstructed files are explicitly named as reconstructed to distinguish recovery material from preserved original inputs.
