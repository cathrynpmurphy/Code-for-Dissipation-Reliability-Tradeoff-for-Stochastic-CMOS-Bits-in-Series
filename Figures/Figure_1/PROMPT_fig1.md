# Figure 1 sources: construction notes

`fig1.pdf` is a pure-vector TikZ/pgfplots rebuild of an earlier
PowerPoint-derived rendering of this figure. The original had three
problems that motivated the rebuild: it embedded 22 raster objects,
including the density panel as a 370x370 indexed image at 218 ppi (below
print quality, with even the math labels rendered as small bitmaps); its
page box held a superseded draft of the figure sitting outside the visible
crop region, so tools reading the page box directly would surface the wrong
version; and no editable source existed for it.

`fig1.tex` is a `figure*`, authored at exactly `\figcanvas{\figfullwidth}{128pt}`
(510pt wide, matching the original figure's 610.4:153.0 aspect ratio at that
width -- `pdfinfo` reports this as 508.09 bp, big points rather than LaTeX
pt). It uses the shared style in `../figstyle.tex`.

## What the figure shows (left to right)

A single wide row, included in the paper at `width=6.5in`:

1. **Left panel, light gray background.** One NOT gate: a PMOS transistor
   above an NMOS transistor, drain rail labelled `$V_{\rm dd}$` (up arrow),
   source rail `$-V_{\rm dd}$` (down arrow), input node `$v_0$` on the left
   as a filled gray-green dot, output node `$v_1$` on the right as a filled
   blue dot. Two pairs of vertical double-headed rate arrows, one pair
   beside the PMOS labelled `$\lambda^{{\rm d},0}_{+}$` /
   `$\lambda^{{\rm u},0}_{+}$` and one beside the NMOS labelled
   `$\lambda^{{\rm d},0}_{-}$` / `$\lambda^{{\rm u},0}_{-}$`.
2. **Middle panel.** The CMOS unit drawn abstractly: two logic-gate
   triangles between the same two nodes `$v_0$` and `$v_1$`, the upper one
   on a light blue background labelled `$k=1$`, the lower one on a light
   gray background labelled `$k=0$`. Dotted leader lines connect this panel
   to panels 1 and 3, showing that each triangle expands into a transistor
   pair.
3. **Right panel, light blue background.** The same as panel 1 but for the
   `$k=1$` gate, with the transistor pair mirrored and the rate arrows on
   the right, labelled with superscript `1` instead of `0`.
4. **Rightmost, the data panel.** A viridis density plot of the
   single-unit steady state over `$v_1$` (horizontal) and `$v_0$`
   (vertical), both axes running -3.2 to 3.2, with a white line along the
   anti-diagonal annotated `$v_0 = -v_1$`. To its right, a narrow attached
   panel showing the bimodal distribution along that anti-diagonal as a
   blue curve on a vertical axis.

## Data for the density panel (panel 4)

`Figure 1 Data/`:

* `Figure1AxisTicks.txt` -- 65 values, -3.2 to 3.2 in steps of 0.1. These
  are node voltages in units of $V_T$; the grid is the $\pm 32 v_e$ cutoff
  with $v_e/V_T = 0.1$.
* `Figure1HeatMap.txt` -- a 65x65 tab-separated matrix, the steady-state
  probability. A handful of entries are tiny negatives (~1e-18), a solver
  artifact clamped to zero during preparation.

`prep_fig1.py` converts these into `fig1_heatmap.dat` (x y z triples for
pgfplots) and `fig1_marginal.dat` (the anti-diagonal cut). The directory
name contains spaces, so it is joined in Python and never handed to LaTeX
directly.

## The density panel is drawn as vector cells

The 65x65 grid is drawn as vector cells (`pgfplots` `matrix plot*` with
`shader=flat corner`), not as an included bitmap -- 4225 filled cells is a
modest file size and compile time. The result is crisper than the original
218-ppi render, though visibly cell-based rather than smoothly
interpolated; that is the deliberate tradeoff of going fully vector.
