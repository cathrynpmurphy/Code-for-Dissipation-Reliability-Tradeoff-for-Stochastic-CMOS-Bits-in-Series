#!/usr/bin/env python3
"""Emit fig5_tau.dat, fig5_dissipation.dat, fig5_tradeoff.dat and fig5_ref.dat.

See FigureSources/README.md. Deterministic; no styling decisions here.

Port of ../PostReviewFigures/make_composite_fig3.py, which is left in place as
the reference. That script carried the whole of the paper's main result as two
literal dicts, DATA_tau and DATA_Q, in its own source. They have been
transcribed verbatim into

    ../PostReviewFigures/Figure 5 Data/tau_err.txt      <tau_err>, in tau_0
    ../PostReviewFigures/Figure 5 Data/dissipation.txt  Qdot, in k_B T / tau_0

(see those files' headers for provenance) and are read from there, so a Zenodo
deposit of this directory plus the raw data reproduces the figure. prep_fig4.py
reads the same tau_err.txt. The directory name has a space in it; it is joined
in Python and never handed to LaTeX.

Outputs:
    fig5_tau.dat          L t11 t12 t13 t14   -- panel (a)
    fig5_dissipation.dat  L q11 q12 q13 q14   -- panel (b)
    fig5_tradeoff.dat     q11 t11 ... q14 t14 -- panel (c), one row per L
    fig5_ref.dat          L y                 -- the two endpoints of the
                          reaction-diffusion reference line of panel (a)
"""
import os
import numpy as np

# See Figures/"Figure 5 Data"/README.txt for what belongs here and why.
DATA = os.path.join("..", "Figure 5 Data")
LS = [1, 3, 5, 7]
VDD = [1.1, 1.2, 1.3, 1.4]


def read_table(path):
    head = [l for l in open(path) if l.lstrip().startswith("#")]
    cols = np.array([float(t) for t in head[-1].split()[2:]])
    tab = np.loadtxt(path)
    return tab[:, 0].astype(int), cols, tab[:, 1:]


L, v1, tau = read_table(os.path.join(DATA, "tau_err.txt"))
L2, v2, Q = read_table(os.path.join(DATA, "dissipation.txt"))
assert list(L) == list(L2) == LS, (L, L2)
assert np.allclose(v1, VDD) and np.allclose(v2, VDD), (v1, v2)
assert tau.shape == Q.shape == (4, 4)
# both quantities rise with L and with V_dd; that is the whole content of (a),(b)
assert (np.diff(tau, axis=0) > 0).all() and (np.diff(tau, axis=1) > 0).all()
assert (np.diff(Q, axis=0) > 0).all() and (np.diff(Q, axis=1) > 0).all()

hdr = lambda p: " ".join(f"{p}{str(v).replace('.','')}" for v in VDD)

np.savetxt("fig5_tau.dat", np.column_stack([L, tau]),
           header="L " + hdr("t"), comments="", fmt=["%d"] + ["%.6g"] * 4)
np.savetxt("fig5_dissipation.dat", np.column_stack([L, Q]),
           header="L " + hdr("q"), comments="", fmt=["%d"] + ["%.6g"] * 4)

# panel (c) pairs the two: one x,y column pair per V_dd, one row per L
cols, names = [], []
for j, v in enumerate(VDD):
    cols += [Q[:, j], tau[:, j]]
    names += [f"q{str(v).replace('.','')}", f"t{str(v).replace('.','')}"]
np.savetxt("fig5_tradeoff.dat", np.column_stack(cols),
           header=" ".join(names), comments="", fmt="%.6g")

# ---- the reaction-diffusion reference line of panel (a) ------------------
# make_composite_fig3.py takes the slope of ln<tau_err> between the first two
# V_dd=1.4 points and continues it straight across the whole L range, so that
# it diverges above the saturating data. It is a straight line on the log
# axis, so two endpoints carry it exactly.
j14 = VDD.index(1.4)
s = (np.log(tau[1, j14]) - np.log(tau[0, j14])) / (LS[1] - LS[0])
ends = np.array([1.0, 7.0])
np.savetxt("fig5_ref.dat", np.column_stack([ends, tau[0, j14] * np.exp(s * (ends - 1))]),
           header="L y", comments="", fmt="%.8g")

print("wrote fig5_tau.dat, fig5_dissipation.dat, fig5_tradeoff.dat, fig5_ref.dat")
print(f"  panel (a) reference slope dln<tau>/dL = {s:.6f} (from V_dd=1.4, L=1->3),"
      f" reaching {tau[0, j14] * np.exp(s * 6):.6g} at L=7")

# matplotlib's autoscale on panel (c), for the limits hard-coded in fig5.tex:
# 5% margins, taken in log space on the log axis.
xr = Q.min(), Q.max()
pad = 0.05 * (xr[1] - xr[0])
lo, hi = np.log10(tau.min()), np.log10(tau.max())
lpad = 0.05 * (hi - lo)
print(f"  panel (c) autoscale: xmin={xr[0]-pad:.4f} xmax={xr[1]+pad:.4f} "
      f"ymin={10**(lo-lpad):.6g} ymax={10**(hi+lpad):.6g}")
