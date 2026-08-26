# Script status and reproducibility notes

The repository preserves the original manuscript-associated Supplementary scripts and, where required, adds narrowly scoped reproducibility wrappers. Original statistical implementations are retained unchanged whenever possible.

## Module status

- `01_alignment_reltime_arpip_iel`: manuscript-associated reference scripts; public input/output/QC files document the available reruns and derived analyses.
- `02_feature_catalogue_registry_qc`: manuscript-associated catalogue and registry validation scripts.
- `03_boundary_pair_universe`: final BPU builder/permutation code is preserved unchanged, but the two historical normalized BPU input CSV files were lost before public archiving. Their recorded SHA-256 checksums and preserved QC/result summaries are documented. The historical BPU is therefore not claimed to be fully rerunnable from the public package.
- `04_pat_aem_coverage_distance`: use the v2 PAT/AEM compatibility wrappers with the original PAT/AEM statistical scripts. Coverage/Distance uses the final manuscript-associated Table S9 runner.
- `05_multiloop_context_dispersion`: use `Ras85D_run_multiloop_dispersion_v2.py` for exact reproduction of the primary fixed-seed analysis; the original v1 file is retained for provenance.

## Historical versus public-rerun files

`SCRIPT_SHA256.csv` records hashes for the original scripts copied from the Supplementary archive. `SCRIPT_PATCH_SHA256.csv` records hashes for the v2 reproducibility files added during public-archive validation. The latter files document compatibility or deterministic-randomization corrections and do not introduce new biological hypotheses or change the reported statistical model.
