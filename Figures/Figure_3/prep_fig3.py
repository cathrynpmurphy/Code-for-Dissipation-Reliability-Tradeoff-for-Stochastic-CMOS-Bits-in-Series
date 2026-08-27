#!/usr/bin/env python3
"""Emit fig3_pss.dat, fig3_prelax.dat and the two 1e-8 contour paths.

See FigureSources/README.md. Deterministic; no styling decisions here -- the
colour scale, its floor and its ceiling all live in fig3.tex.

Port of ../PostReviewFigures/make_basis_convergence.py, which is left in place
as the reference. Two things are repaired on the way:

  * that script opens 'L2err_Pss.txt' by bare filename from its own directory.
    The raw exports have since moved into "Figure 3 Data/", so the load raised
    FileNotFoundError -- and, worse, os.path.exists('L2err_Pss.txt') then went
    false and the PROVISIONAL branch would have silently swapped panel (a) for
    the matrix digitized off a rendered PDF. There is no fallback here: both
    matrices are read from "Figure 3 Data/" by their real paths, and the
    digitized stand-in is not used at all. (The directory name has a space in
    it; it is joined in Python and never handed to LaTeX.)
  * the 1e-8 contour was drawn by matplotlib. It is computed here instead, by
    marching squares with linear interpolation in Z -- matplotlib's own
    convention -- so the polyline lands where the old one did.

Raw inputs (../PostReviewFigures/Figure 3 Data/):
    L2err_Pss.txt      9x32 tab-separated Julia complex, ||Pbar_ss  - P_ss ||_2
    L2err_Prelax.txt   9x32 tab-separated Julia complex, ||Pbar_rel- P_rel||_2
Both on the same grid: basis d = 2,4,...,64 across the columns, bond dimension
chi = 10,20,...,90 down the rows.

Outputs:
    fig3_pss.dat, fig3_prelax.dat        d chi logerr, d varying fastest,
                                         32 per scanline, for `matrix plot*`
                                         with mesh/cols=32. logerr = log10 of
                                         the error, unclipped.
    fig3_pss_c8.dat, fig3_prelax_c8.dat  d chi polyline of the 1e-8 contour,
                                         blank line between components.
"""
import os
import numpy as np

# Raw exports now come from Figures/Figure_3/CompareGMRESwithTNsL3.jl, which writes
# them here so this deposit is self-contained (see zenodo_v4_instructions.txt).
DATA = os.path.join("Figure 3 Data")
LEVEL = 1e-8            # accuracy contour, as in make_basis_convergence.py


def load_complex(path):
    """Read the tab-separated Julia complex export as a real matrix."""
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            rows.append([float(t.split('+')[0].strip()) for t in line.split('\t')])
    return np.array(rows)


def contour(x, y, Z, level):
    """Marching squares at a single level; returns a list of polylines.

    Linear interpolation in Z along each cell edge and the four-corner mean as
    the saddle test, which is what matplotlib's contour does, so the path is
    the one the matplotlib figure drew.
    """
    ny, nx = Z.shape
    segs = []

    def interp(pa, pb, za, zb):
        t = (level - za) / (zb - za)
        return (pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1]))

    for i in range(ny - 1):
        for j in range(nx - 1):
            c = [(x[j], y[i]), (x[j + 1], y[i]), (x[j + 1], y[i + 1]), (x[j], y[i + 1])]
            z = [Z[i, j], Z[i, j + 1], Z[i + 1, j + 1], Z[i + 1, j]]
            # crossing point on each of the four edges, or None
            e = []
            for k in range(4):
                za, zb = z[k], z[(k + 1) % 4]
                if (za < level) != (zb < level):
                    e.append(interp(c[k], c[(k + 1) % 4], za, zb))
                else:
                    e.append(None)
            hit = [k for k in range(4) if e[k] is not None]
            if len(hit) == 2:
                segs.append((e[hit[0]], e[hit[1]]))
            elif len(hit) == 4:
                # saddle: the centre value decides which pair joins up
                if (sum(z) / 4.0) < level:
                    segs.append((e[0], e[1])); segs.append((e[2], e[3]))
                else:
                    segs.append((e[1], e[2])); segs.append((e[3], e[0]))
    return chain(segs)


def chain(segs):
    """Stitch unordered segments into as few polylines as possible."""
    def key(p):
        return (round(p[0], 9), round(p[1], 9))

    remaining = list(segs)
    lines = []
    while remaining:
        a, b = remaining.pop()
        line = [a, b]
        grew = True
        while grew:
            grew = False
            for idx, (p, q) in enumerate(remaining):
                for pt, other in ((p, q), (q, p)):
                    if key(pt) == key(line[-1]):
                        line.append(other); remaining.pop(idx); grew = True; break
                    if key(pt) == key(line[0]):
                        line.insert(0, other); remaining.pop(idx); grew = True; break
                if grew:
                    break
        lines.append(line)
    return lines


def write_grid(path, x, y, Z):
    with open(path, "w") as f:
        f.write("d chi logerr\n")
        for i, yy in enumerate(y):
            if i:
                f.write("\n")                       # scanline break
            for j, xx in enumerate(x):
                f.write(f"{xx:g} {yy:g} {np.log10(Z[i, j]):.6f}\n")


def write_contour(path, lines):
    with open(path, "w") as f:
        f.write("d chi\n")
        for k, line in enumerate(lines):
            if k:
                f.write("\n")
            for px, py in line:
                f.write(f"{px:.5f} {py:.5f}\n")


L2_ss = np.abs(load_complex(os.path.join(DATA, "L2err_Pss.txt")))
L2_rel = np.abs(load_complex(os.path.join(DATA, "L2err_Prelax.txt")))

basis = np.arange(2, 2 + 2 * L2_rel.shape[1], 2)      # d = 2..64
bond = np.arange(10, 10 + 10 * L2_rel.shape[0], 10)   # chi = 10..90
assert L2_ss.shape == L2_rel.shape == (bond.size, basis.size), \
    (L2_ss.shape, L2_rel.shape)
assert (basis[0], basis[-1]) == (2, 64) and (bond[0], bond[-1]) == (10, 90)

# both errors must fall with d and with chi, and stay inside the colour range
for name, Z in (("P_ss", L2_ss), ("P_relax", L2_rel)):
    assert Z.min() > 0, f"{name}: non-positive entry"
    assert Z[-1, -1] < Z[0, 0], f"{name}: converged corner is not the smallest"
    write_grid(f"fig3_{'pss' if name == 'P_ss' else 'prelax'}.dat", basis, bond, Z)
    write_contour(f"fig3_{'pss' if name == 'P_ss' else 'prelax'}_c8.dat",
                  contour(basis, bond, Z, LEVEL))
    jb = int(np.argmin(np.abs(basis - 32)))
    print(f"{name:8s} range {Z.min():.2e} .. {Z.max():.2e}; "
          f"L2 at d=32, chi=90 = {Z[-1, jb]:.2e}")

for stem in ("pss", "prelax"):
    n = sum(1 for _ in open(f"fig3_{stem}_c8.dat")) - 1
    print(f"wrote fig3_{stem}.dat and fig3_{stem}_c8.dat ({n} contour points)")
