# Figures 3, 4 and 5 sources: construction notes

`fig3.pdf`, `fig4.pdf` and `fig5.pdf` are pure-vector TikZ/pgfplots ports of
figures originally produced in matplotlib. Unlike Figures 1 and 2, these
three already had working generator scripts (`make_basis_convergence.py`,
`make_fig4_vdd2.py`, `make_composite_fig3.py`, kept alongside as the
reference); the port preserves every substantive choice -- panel
arrangement, what is plotted against what, colour meaning, contour levels,
annotations -- and changes only the rendering layer. All three use the
shared style in `../figstyle.tex`.

## Two corrections made during the port

1. **Data paths.** `make_basis_convergence.py` originally loaded
   `L2err_Pss.txt`, `L2err_Prelax.txt` and `Pss_L2_digitized.txt` by bare
   filename from its own directory; those files now live in
   `Figure_3/"Figure 3 Data"/`, alongside the ported script, so `prep_fig3.py`
   reads them from there. The directory name contains spaces, so it is
   joined in Python and never handed to LaTeX directly.
2. **Production numbers moved out of source.** `make_composite_fig3.py`
   originally carried the results of the paper's main calculation as
   literal dicts in its source: `DATA_tau` (mean time-to-error, units of
   $\tau_0$) and `DATA_Q` (steady-state dissipation, units of
   $k_{\rm B}T/\tau_0$), keyed by $V_{\rm dd}$ over $L = 1,3,5,7$. For this
   deposit, that is the wrong place for them: they are extracted verbatim
   into `Figures/"Figure 5 Data"/tau_err.txt` and `dissipation.txt`, with a
   header recording that they were transcribed from
   `make_composite_fig3.py` and originate from the paper's production DMRG
   runs. `prep_fig4.py` and `prep_fig5.py` both read from these two files,
   so Figures 4 and 5 cannot disagree about the production numbers.

## Figure 3 -- basis and bond-dimension convergence

Source: `make_basis_convergence.py`. Two stacked panels at $L=3$ on a
shared grid and a single shared log colour scale with one colorbar: (a)
$P_{\rm ss}$, (b) $P_{\rm relax}$, both the $L_2$ eigenvector error. Axes:
number of retained basis vectors $d$ horizontal, bond dimension $\chi$
vertical. Grid is $d = 2,4,\ldots,64$ and $\chi = 10,20,\ldots,90$; the
matrices are 9x32, $\chi$ down the rows and $d$ across the columns. The
white $10^{-8}$ accuracy contour and the red dashed production cutoff at
$d=32$ are drawn explicitly, with white-boxed panel letters.

Both density panels and the colorbar were originally rasters at 100 ppi,
the lowest-resolution art in the paper; the grid is small (9x32), so
drawing the cells as vector fills carries no size cost.

Single column, `\figcanvas{\figcolwidth}{<h>}` with height around 250pt,
included with `[width=\columnwidth]`.

## Figure 4 -- Arrhenius scaling in $V_{\rm dd}^2$

Source: `make_fig4_vdd2.py`. Its data comes from the same production
numbers as Figure 5 (see above), so both read from the same extracted
files rather than duplicating them.

$\ln\langle\tau_{\rm err}\rangle$ against $V_{\rm dd}^2$, one series per
chain length $L = 1,3,5,7$ with markers and linear fits. The paper quotes
$R^2 \geq 0.997$ for every chain; the fits are recomputed from the
extracted data in the port, and the actual values are reported by
`prep_fig4.py`.

Single column, height around 190pt, `[width=\columnwidth]`.

## Figure 5 -- the main results figure

Source: `make_composite_fig3.py`. Three panels: (a) mean bit-flip time
against $L$ on a log axis, one curve per $V_{\rm dd}$, with a dashed
reference line showing what pure exponential-in-$L$ growth would look
like; (b) dissipation $\dot{Q}$ against $L$; (c) the tradeoff,
$\langle\tau_{\rm err}\rangle$ against $\dot{Q}$. The $V_{\rm dd}$ colour
ramp means the same thing in all three panels and matches Figure 4's.

Full width, `\figcanvas{\figfullwidth}{<h>}` with height around 200pt,
included with `[width=\textwidth]`.
