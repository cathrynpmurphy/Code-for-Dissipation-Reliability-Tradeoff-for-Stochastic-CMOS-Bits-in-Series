#!/usr/bin/env python3
"""Emit fig2_basis.dat from the raw exports (see FigureSources/README.md)."""
import os
import numpy as np

# Raw exports now come from GetEigenvectorsWithDMRG/ExcitedStateDMRG.jl, which
# writes them here so this deposit is self-contained (see zenodo_v4_instructions.txt).
DATA = os.path.join("Figure 2 Data")

x = np.loadtxt(os.path.join(DATA, "BasisVectorsXAxis.txt"))
y = np.loadtxt(os.path.join(DATA, "BasisVectorsYAxis.txt"))
assert x.shape == (61,), x.shape
assert y.shape == (5, 61), y.shape

# parity ladder check: row j should have j sign changes, parity (-1)^j
for j, row in enumerate(y):
    s = np.sign(row)[np.sign(row) != 0]
    assert int((np.diff(s) != 0).sum()) == j, f"row {j}: wrong node count"
    assert np.sign(np.corrcoef(row, row[::-1])[0, 1]) == (-1) ** j, f"row {j}: wrong parity"

np.savetxt("fig2_basis.dat", np.column_stack([x, *y]),
           header="x f0 f1 f2 f3 f4", comments="", fmt="%.8g")
print("wrote fig2_basis.dat")
