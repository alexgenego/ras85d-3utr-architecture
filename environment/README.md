# Computational environment

This directory documents the software environment used to validate the archived custom Python workflows associated with the manuscript.

## Python environment

The archived scripts were validation-tested with:

- Python 3.13.5
- NumPy 2.3.5
- pandas 2.2.3
- SciPy 1.17.0
- statsmodels 0.14.6
- openpyxl 3.1.5
- NetworkX 3.6.1

`requirements.txt` records these package versions. `environment_check.py` can be used to print the Python and package versions available in a local environment.

`openpyxl` is required as the Excel engine used by pandas for manuscript-associated `.xlsx` inputs. `statsmodels` is required by the branch-statistics workflow. NetworkX is retained in the archived analytical environment because it was used during network-oriented processing and validation, although not every archived script imports it directly.

## External software and resources

Additional stages of the study used external software and databases described in the manuscript and Supplementary Materials, including MEGA, ARPIP/Bio++, MEME/MAST, CENSOR/Repbase, TargetScanFly, and RNAfold/ViennaRNA. These programs are not installed by `requirements.txt`.

Fixed random seeds, permutation counts, analysis-specific parameters, and quality-control information are documented in the corresponding scripts, README files, and `qc/` records.
