#!/usr/bin/env python3
"""Emit fig1_heatmap.dat and fig1_marginal.dat from the raw exports.

See FigureSources/README.md. Deterministic; no styling decisions here.

Raw inputs (directory name contains spaces -- joined in Python, never handed
to LaTeX):
    Figure1AxisTicks.txt  65 node voltages, -3.2 .. 3.2 step 0.1, units of V_T.
                          The grid is the +-32 v_e cutoff with v_e/V_T = 0.1.
    Figure1HeatMap.txt    65x65 tab-separated steady-state probability.

Outputs:
    fig1_heatmap.dat   x y z triples, x (=v1) varying fastest, 65 per scanline,
                       for pgfplots `matrix plot*` with mesh/cols=65.
    fig1_marginal.dat  the cut along the anti-diagonal v0 = -v1, as v0 vs p,
                       p normalised to unit peak.

Index convention: M[i][j] = P(v0 = t[i], v1 = t[j]) -- the first (row) index is
v0, which the figure puts on the vertical axis. The matrix is symmetric to
within 4.0e-5, i.e. 0.6% of its peak, and that residue sits entirely in the
dark tails, so the transpose is indistinguishable in the rendered panel; the
original 218-ppi bitmap is interpolated and cannot settle it either. Set
TRANSPOSE = True to flip if the solver's convention is ever confirmed to be
the other way round.
"""
import os
import numpy as np

TRANSPOSE = False

# Raw exports now come from GetL3EigenvectorsWithGMRES/SparseSolve.jl, which
# writes them here so this deposit is self-contained (see zenodo_v4_instructions.txt).
DATA = os.path.join("Figure 1 Data")

t = np.loadtxt(os.path.join(DATA, "Figure1AxisTicks.txt"))
M = np.loadtxt(os.path.join(DATA, "Figure1HeatMap.txt"))
assert t.shape == (65,), t.shape
assert M.shape == (65, 65), M.shape

# the grid is the +-32 v_e cutoff, v_e/V_T = 0.1
assert np.allclose(t, np.arange(-32, 33) * 0.1), "axis ticks are not -3.2..3.2 step 0.1"

# the solver leaks ~1e-18 negatives; clamp at zero, and check that is all it is
assert M.min() > -1e-12, f"unexpectedly large negative entry {M.min():g}"
assert abs(M.sum() - 1.0) < 1e-9, f"not normalised: sum = {M.sum():.12g}"
M = np.clip(M, 0.0, None)

if TRANSPOSE:
    M = M.T

# the distribution is supported on the anti-diagonal: the two modes must sit at
# (v0, v1) = (-+1.1, +-1.1), i.e. mirror images through the origin
i, j = np.unravel_index(np.argmax(M), M.shape)
assert (t[i], t[j]) == (1.1, -1.1) or (t[i], t[j]) == (-1.1, 1.1), (t[i], t[j])
assert np.isclose(M[i, j], M[64 - i, 64 - j]), "modes are not symmetric through the origin"

# ---- density panel ------------------------------------------------------
# x = v1 (horizontal), y = v0 (vertical), x varying fastest
rows = [(t[jj], t[ii], M[ii, jj]) for ii in range(65) for jj in range(65)]
with open("fig1_heatmap.dat", "w") as f:
    f.write("x y z\n")
    for k, (x, y, z) in enumerate(rows):
        if k and k % 65 == 0:
            f.write("\n")            # scanline break, 65 columns per row
        f.write(f"{x:.4g} {y:.4g} {z:.10g}\n")
print(f"wrote fig1_heatmap.dat ({len(rows)} cells, peak {M.max():.6g})")

# ---- marginal along the anti-diagonal v0 = -v1 --------------------------
# v1 = -v0 is exactly on-grid because t is symmetric: t[64-k] == -t[k]
assert np.allclose(t[::-1], -t)
cut = np.array([M[k, 64 - k] for k in range(65)])
with open("fig1_marginal.dat", "w") as f:
    f.write("v0 p\n")
    for v, p in zip(t, cut / cut.max()):
        f.write(f"{v:.4g} {p:.8g}\n")
print(f"wrote fig1_marginal.dat (bimodal peaks at v0 = "
      f"{t[cut.argmax()]:+.1f} and {-t[cut.argmax()]:+.1f})")
