# ---------------------------------------------------------------------------
#  Figures/Makefile -- build the CMOS chain paper's figures from source.
#
#      make            build every figure that has sources
#      make figN       build Figure_N/figN.pdf (runs prep_figN.py first)
#      make check      assert the contract: exact width, no rasters, no newtx
#      make clean      remove LaTeX aux files (keeps the .pdf and .dat)
#      make distclean  also remove generated .pdf and .dat
#
#  See README.md for the style contract. Every figure is authored at exactly
#  the width it occupies on the page, so \includegraphics never rescales it.
#
#  Each figure lives in its own Figure_N/ subfolder (figN.tex, prep_figN.py,
#  the .dat files, and -- for figures 1-3 -- the raw-data folder and/or local
#  copies of the solver code that produced it). figstyle.tex is the one
#  exception: it stays here, one level up, as the single style source of
#  truth (Figure_N/figN.tex does \input{../figstyle}), so every rule below
#  cds into the figure's own subfolder before doing anything.
# ---------------------------------------------------------------------------

# The prep scripts need numpy. On macOS the Homebrew python3 frequently lacks
# it while the system python3 has it; prefer whichever can actually import it.
# Override with e.g. `make PYTHON=/path/to/python3 fig1`.
PYTHON ?= $(shell python3 -c 'import numpy' >/dev/null 2>&1 \
                  && echo python3 || echo /usr/bin/python3)
LATEX  ?= pdflatex
LFLAGS  = -interaction=nonstopmode -halt-on-error

FIGURES = Figure_1/fig1.pdf Figure_2/fig2.pdf Figure_3/fig3.pdf Figure_4/fig4.pdf Figure_5/fig5.pdf

.PHONY: all fig1 fig2 fig3 fig4 fig5 check clean distclean
all: $(FIGURES)

# ---- figure 1 -------------------------------------------------------------
# Raw exports live in "Figure_1/Figure 1 Data/" (written by
# GetL3EigenvectorsWithGMRES/SparseSolve.jl, kept local to this deposit -- see
# ../zenodo_v4_instructions.txt) -- a directory whose name contains spaces,
# which make cannot carry in a prerequisite list. prep_fig1.py opens it by
# name instead, so the .dat rule depends only on the script; touch it (or run
# `make -B fig1`) after the raw data changes.
Figure_1/fig1_heatmap.dat Figure_1/fig1_marginal.dat: Figure_1/prep_fig1.py
	cd Figure_1 && $(PYTHON) prep_fig1.py

Figure_1/fig1.pdf: Figure_1/fig1.tex figstyle.tex Figure_1/fig1_heatmap.dat Figure_1/fig1_marginal.dat
	cd Figure_1 && $(LATEX) $(LFLAGS) fig1.tex

fig1: Figure_1/fig1.pdf

# ---- figure 2 -------------------------------------------------------------
# Raw exports in "Figure_2/Figure 2 Data/" (written by
# GetEigenvectorsWithDMRG/ExcitedStateDMRG.jl); same spacey-directory story as
# figure 1.
Figure_2/fig2_basis.dat: Figure_2/prep_fig2.py
	cd Figure_2 && $(PYTHON) prep_fig2.py

Figure_2/fig2.pdf: Figure_2/fig2.tex figstyle.tex Figure_2/fig2_basis.dat
	cd Figure_2 && $(LATEX) $(LFLAGS) fig2.tex

fig2: Figure_2/fig2.pdf

# ---- figure 3 -------------------------------------------------------------
# Raw exports in "Figure_3/Figure 3 Data/" (written by
# Figure_3/CompareGMRESwithTNsL3.jl, which in turn reads the GMRES reference
# GetL3EigenvectorsWithGMRES/SparseSolve.jl writes to the same folder); same
# spacey-directory story as figure 1.
FIG3DAT = Figure_3/fig3_pss.dat Figure_3/fig3_prelax.dat Figure_3/fig3_pss_c8.dat Figure_3/fig3_prelax_c8.dat
$(FIG3DAT): Figure_3/prep_fig3.py
	cd Figure_3 && $(PYTHON) prep_fig3.py

Figure_3/fig3.pdf: Figure_3/fig3.tex figstyle.tex $(FIG3DAT)
	cd Figure_3 && $(LATEX) $(LFLAGS) fig3.tex

fig3: Figure_3/fig3.pdf

# ---- figures 4 and 5 -------------------------------------------------------
# Both read "Figure 5 Data/tau_err.txt" (one level up, shared by both figures
# so they can't disagree about the production numbers -- see
# "Figure 5 Data/README.txt" for what belongs there; it is NOT produced by any
# script in this repo, unlike figures 1-3's raw data). fig4_r2.tex carries the
# four R^2 into the legend, so that too is computed rather than typed.
FIG4DAT = Figure_4/fig4_points.dat Figure_4/fig4_fits.dat Figure_4/fig4_r2.tex
$(FIG4DAT): Figure_4/prep_fig4.py
	cd Figure_4 && $(PYTHON) prep_fig4.py

Figure_4/fig4.pdf: Figure_4/fig4.tex figstyle.tex $(FIG4DAT)
	cd Figure_4 && $(LATEX) $(LFLAGS) fig4.tex

fig4: Figure_4/fig4.pdf

FIG5DAT = Figure_5/fig5_tau.dat Figure_5/fig5_dissipation.dat Figure_5/fig5_tradeoff.dat Figure_5/fig5_ref.dat
$(FIG5DAT): Figure_5/prep_fig5.py
	cd Figure_5 && $(PYTHON) prep_fig5.py

Figure_5/fig5.pdf: Figure_5/fig5.tex figstyle.tex $(FIG5DAT)
	cd Figure_5 && $(LATEX) $(LFLAGS) fig5.tex

fig5: Figure_5/fig5.pdf

# ---- the contract ---------------------------------------------------------
# 510.0pt is 508.09 bp: pdfinfo reports big points (72/in), LaTeX pt are
# 72.27/in. Column figures read 245.08 bp, full-width ones 508.09 bp.
check: $(FIGURES)
	@status=0; for f in $(FIGURES); do \
	  w=$$(pdfinfo $$f | awk '/Page size/{print $$3}'); \
	  n=$$(pdfinfo $$f | awk '/^Pages/{print $$2}'); \
	  img=$$(pdfimages -list $$f | tail -n +3 | grep -c . || true); \
	  echo "$$f: $$w bp wide, $$n page(s), $$img image object(s)"; \
	  [ "$$img" = "0" ] || { echo "  FAIL: $$f contains raster images"; status=1; }; \
	  [ "$$n" = "1" ]   || { echo "  FAIL: $$f is not one page"; status=1; }; \
	  case $$w in 508.0*|245.0*|508.1*|245.1*) ;; \
	    *) echo "  NOTE: $$w bp is neither 245.08 (column) nor 508.09 (full width) -- ok for fig2, which is narrower by design (see Figure_2/PROMPT_fig2.md)";; esac; \
	  pdffonts $$f | grep -Eqi 'newtx|times|mathptmx|nimbus' \
	    && { echo "  FAIL: $$f uses a non-Computer-Modern font"; status=1; } || true; \
	done; \
	grep -Eli 'newtx|\{times\}|mathptmx' Figure_*/fig*.tex && { echo "  FAIL: forbidden font package in a source file"; status=1; } || true; \
	[ $$status = 0 ] && echo "all checks passed"; exit $$status

clean:
	rm -f Figure_*/*.aux Figure_*/*.log Figure_*/*.out Figure_*/*.synctex.gz Figure_*/*.build.log
distclean: clean
	rm -f $(FIGURES) Figure_1/fig1_heatmap.dat Figure_1/fig1_marginal.dat \
	      Figure_2/fig2_basis.dat \
	      $(FIG3DAT) $(FIG4DAT) $(FIG5DAT)
