#!/usr/bin/env python3
"""Emit fig4_points.dat, fig4_fits.dat and fig4_r2.tex.

See FigureSources/README.md. Deterministic; no styling decisions here.

Port of ../PostReviewFigures/make_fig4_vdd2.py, which is left in place as the
reference. That script carried the production time-to-error matrix as a literal
dict in its source; the numbers now live in
../PostReviewFigures/Figure 5 Data/tau_err.txt (transcribed verbatim, see the
header of that file) and are read from there, so Fig. 4 and Fig. 5 cannot drift
apart. The directory name has a space in it; it is joined in Python and never
handed to LaTeX.

The fits are ordinary least squares of ln<tau_err> on V_dd^2, one per chain
length, exactly as in the reference script (np.polyfit degree 1). The paper
claims R^2 >= 0.997 for every chain; that is asserted here, and the values are
written to fig4_r2.tex so the legend cannot fall out of step with the data.

Outputs:
    fig4_points.dat  v2 y1 y3 y5 y7   -- ln<tau_err> at the four V_dd^2
    fig4_fits.dat    v2 f1 f3 f5 f7   -- the fitted lines, two endpoints each
    fig4_r2.tex      \def\RsqOne{...} etc., \input by fig4.tex
"""
import os
import numpy as np

# Shared with prep_fig5.py -- both figures read the same transcribed production
# numbers, kept in one place directly under Figures/ so they cannot drift apart
# (see zenodo_v4_instructions.txt and Figures/"Figure 5 Data"/README.txt).
DATA = os.path.join("..", "Figure 5 Data")
LS = [1, 3, 5, 7]
PAD = 0.03            # the fit lines overhang the data by this much in V_dd^2,
                      # as in make_fig4_vdd2.py
WORDS = {1: "One", 3: "Three", 5: "Five", 7: "Seven"}


def read_table(path):
    """Return (L column, V_dd column labels, values) from one of the production tables."""
    head = [l for l in open(path) if l.lstrip().startswith("#")]
    cols = np.array([float(t) for t in head[-1].split()[2:]])
    tab = np.loadtxt(path)
    return tab[:, 0].astype(int), cols, tab[:, 1:]


L, Vdd, tau = read_table(os.path.join(DATA, "tau_err.txt"))
assert list(L) == LS, L
assert np.allclose(Vdd, [1.1, 1.2, 1.3, 1.4]), Vdd
assert tau.shape == (4, 4) and (tau > 0).all()

V2 = Vdd ** 2
lntau = np.log(tau)

fits, r2 = {}, {}
for i, l in enumerate(LS):
    y = lntau[i]
    m, b = np.polyfit(V2, y, 1)
    yh = m * V2 + b
    fits[l] = (m, b)
    r2[l] = 1 - np.sum((y - yh) ** 2) / np.sum((y - y.mean()) ** 2)
    assert m > 0, f"L={l}: fitted slope is not positive"

# the paper's claim, in Sec. V A: "R^2 >= 0.997 for every chain length studied"
worst = min(r2.values())
assert worst >= 0.997, f"paper claims R^2 >= 0.997; worst fit is {worst:.6f}"

xf = np.array([V2[0] - PAD, V2[-1] + PAD])

with open("fig4_points.dat", "w") as f:
    f.write("v2 " + " ".join(f"y{l}" for l in LS) + "\n")
    for j, x in enumerate(V2):
        f.write(f"{x:.6g} " + " ".join(f"{lntau[i, j]:.8g}" for i in range(4)) + "\n")

with open("fig4_fits.dat", "w") as f:
    f.write("v2 " + " ".join(f"f{l}" for l in LS) + "\n")
    for x in xf:
        f.write(f"{x:.6g} " + " ".join(f"{fits[l][0] * x + fits[l][1]:.8g}"
                                      for l in LS) + "\n")

with open("fig4_r2.tex", "w") as f:
    f.write("% written by prep_fig4.py -- do not edit; R^2 of the four fits\n")
    for l in LS:
        f.write(f"\\def\\Rsq{WORDS[l]}{{{r2[l]:.3f}}}\n")

print("wrote fig4_points.dat, fig4_fits.dat, fig4_r2.tex")
for l in LS:
    m, b = fits[l]
    print(f"  L={l}: ln<tau_err> = {m:.4f} V_dd^2 + {b:+.4f},  R^2 = {r2[l]:.6f}")
print(f"  worst R^2 = {worst:.6f}  (paper claims >= 0.997)")
