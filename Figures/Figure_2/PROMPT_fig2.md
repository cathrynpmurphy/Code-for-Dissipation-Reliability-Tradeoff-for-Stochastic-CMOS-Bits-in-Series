# Figure 2 sources: construction notes

`fig2.pdf` is a pure-vector TikZ/pgfplots rebuild of an earlier
PowerPoint-derived rendering of this figure. Unlike Figure 1, that original
was already fully vector, so the rebuild was mainly about fonts, source
control, and four content corrections (below), rather than about
vectorness: its fonts were Cambria Math, Arial, Aptos and
DejaVuSans-ExtraLight, none of them Computer Modern, so no text in the
figure matched the rest of the paper; and its page box held a superseded
draft of the whole diagram sitting outside the visible crop region, with
`$v_1 \ldots v_L$` instead of `$v_0 \ldots v_L$` and
`$\mu_1 = \tau_{err}^{-1}$`, which contradicts
$\langle\tau_{\rm err}\rangle = 2/\mu_1$ in the body -- a discrepancy the
corrections below resolve.

The figure targets the full 246pt column, height 300-380pt, re-proportioned
from the original's more severely vertical layout (204x427pt) so it reads
comfortably at that width; `\figcanvas{<width>}{<height>}` pins the canvas,
and a correctly-sized 246pt-wide figure reads 245.08 bp in `pdfinfo`. It
uses the shared style in `../figstyle.tex` and follows Figure 1's idioms
for circuit symbols and node dots so the two read as siblings.

## What the figure shows (top to bottom, single column, 2.25in wide)

1. A chain of NOT-gate triangle pairs with nodes labelled `$v_0$`, `$v_i$`,
   `$v_{i+1}$`, `$v_L$`, ellipses between them.
2. A row of gray circles (the MPS tensors) sitting above a row of blue
   rounded squares (the MPO tensors), joined by vertical lines, with
   horizontal bond lines and short downward free legs. The rightmost site
   is picked out by a pale blue highlight box.
3. The eigenvalue equation
   `$\hat{\mathbb{W}}|P_i(\mathbf{v})\rangle = \mu_i|P_i(\mathbf{v})\rangle$`
   with a downward arrow.
4. A row of gray circles alone (the resulting MPS), rightmost one darker
   and highlighted.
5. A zoom callout from that highlighted site to an inset plot of the
   single-unit basis vectors (see data below).
6. Two arrows labelled `DMRG` and `Excited State DMRG`, descending from the
   eigenvalue-equation / MPS stage (items 3-4) so the eye reads "DMRG
   produces these states," with the basis-function inset hanging off to
   the side as the callout it is.
7. Two schematic curve panels at the bottom, titled
   `$|P_{\rm ss}\rangle$` and `$|P_{\rm relax}\rangle$` so the arrow labels
   read as the methods rather than as titles for the curves: a symmetric
   double-lobed one marked `$\mu_0 = 0$`, and an antisymmetric one marked
   `$\mu_1 = 2\langle\tau_{\rm err}\rangle^{-1}$` with a red annotation
   `Probability transfer` pointing at the negative lobe. These panels are
   deliberately schematic -- for $L>1$ these eigenvectors are functions of
   all $L+1$ node voltages and cannot be plotted directly, so they are
   drawn as smooth analytic curves with no numeric axes, which would
   otherwise imply they were data.

The inset equation is 0-indexed to match the body text's convention
($\phi_0 = P_{\rm ss}^{(L=1)}$, retained set $\{\varphi_i\}_{i=0}^{d-1}$):
`$c_0\varphi_0 + c_1\varphi_1 + \ldots c_{d-1}\varphi_{d-1}$`, with
$\varphi_0$ and $\varphi_1$ (the two curves drawn in blue) highlighted.

## Data for the inset (item 5)

`Figure 2 Data/`:

* `BasisVectorsXAxis.txt` -- 61 values, -30 to 30, i.e. $v_L/v_e$.
* `BasisVectorsYAxis.txt` -- 5 tab-separated rows of 61 values: the
  marginalized, orthonormalized single-unit basis vectors
  $\varphi_0 \ldots \varphi_4$. Row 0 is symmetric with no sign change (the
  steady state), row 1 antisymmetric with one, and so on up to four sign
  changes -- a clean parity ladder, which also serves as a correctness
  check on the loading code.

`prep_fig2.py` emits `fig2_basis.dat` with columns `x f0 f1 f2 f3 f4`. The
directory name contains spaces, so it is joined in Python and never handed
to LaTeX directly. Style: $\varphi_0$ and $\varphi_1$ as thick blue curves
(a dark and a mid blue), $\varphi_2$--$\varphi_4$ as thin pale green. Axis
labels `$v_L/v_e$` and `$\varphi_j(v_L/v_e)$`.

`fig2_inset_reference.tex` compiles standalone with zero errors, zero
raster objects, and Computer Modern only, and is the starting point folded
into `fig2.tex` for the inset panel; `prep_fig2.py` asserts the parity
ladder as a load check on `fig2_basis.dat`, confirming the rows are loaded
in the right order.
