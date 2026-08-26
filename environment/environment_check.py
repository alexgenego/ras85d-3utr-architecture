#!/usr/bin/env python3
"""Print versions of Python packages used in the archived Ras85D workflows."""

from __future__ import annotations

import platform

import networkx
import numpy
import openpyxl
import pandas
import scipy
import statsmodels


PACKAGES = (
    ("NumPy", numpy.__version__),
    ("pandas", pandas.__version__),
    ("SciPy", scipy.__version__),
    ("statsmodels", statsmodels.__version__),
    ("openpyxl", openpyxl.__version__),
    ("NetworkX", networkx.__version__),
)

print(f"Python {platform.python_version()}")
for name, version in PACKAGES:
    print(f"{name} {version}")
