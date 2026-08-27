# Code for "Dissipation-Reliability Tradeoff for Stochastic CMOS Bits in Series"

This repository contains Julia code for computing the dissipation-reliability
tradeoff in stochastic CMOS circuits using tensor network methods, as
described in the associated publication. The code is written in Julia and
relies primarily on [ITensors.jl](https://itensor.org) (Fishman, White, and
Stoudenmire, *SciPost Phys. Codebases* 2022).

This README shows:: for every
figure in the paper, this document walks from "I want Figure N" back through
the script that drew it, the data it plots, the script that produced that
data, and the software that produced the raw data underneath it.

## Repository layout

| Folder | What it does |
|---|---|
| `Solvers/GetEigenvectorsWithDMRG/` | Tensor-network (DMRG) solver: builds the chain's master equation as an MPO, applies the L=1 basis transformation, and warm-starts DMRG up to the production basis size/bond dimension to get the steady state (`SteadyStateDMRG.jl`) and first-excited/relaxation eigenvector (`ExcitedStateDMRG.jl`). |
| `Solvers/GetL3EigenvectorsWithGMRES/` | Independent, exact/sparse reference solve at L=3: builds the same rate matrix directly as a sparse matrix (not an MPO) and diagonalizes it with `KrylovKit`. This is the ground truth Figure 3 checks the tensor-network method against. |
| `Solvers/ComputeDissipation/` | Given a converged steady-state MPS, builds the dissipation MPO and computes the steady-state dissipation rate $\dot Q$. |
| `Figures/` | One subfolder per figure (`Figure_1` .. `Figure_5`), each containing the TikZ/pgfplots source, the prep script that turns raw data into plotting-ready `.dat` files, and (where applicable) the Julia script and/or local copies of the solver code that produced the raw data. See `Figures/README.md` for the figure style contract. |

## The production sweep and its Data/ folders

The paper's results come from 16 production points: $L = 1,3,5,7$ (number of
chain units) crossed with $V_{\rm dd}/V_T = 1.1,1.2,1.3,1.4$. Both
`Solvers/GetEigenvectorsWithDMRG/SteadyStateDMRG.jl` and
`Solvers/GetEigenvectorsWithDMRG/ExcitedStateDMRG.jl`, plus
`Solvers/ComputeDissipation/ComputeDissipation.jl`, now take `L` and `vDD_vT` as
command-line arguments (`julia SteadyStateDMRG.jl <L> <vDD_vT>`; omitting them
falls back to L=7, vDD_vT=1.2, matching what used to be hard-coded), and each
writes into a `Data/` folder alongside itself, named for the point it was run
at:

```
Solvers/GetEigenvectorsWithDMRG/Data/Pss_Vdd<vDD_vT>_L<L>.h5        (SteadyStateDMRG.jl)
Solvers/GetEigenvectorsWithDMRG/Data/Prelax_Vdd<vDD_vT>_L<L>.h5      (ExcitedStateDMRG.jl)
Solvers/GetEigenvectorsWithDMRG/Data/mu1_Vdd<vDD_vT>_L<L>.jld2       (ExcitedStateDMRG.jl)
Solvers/ComputeDissipation/Data/Dissipation_Vdd<vDD_vT>_L<L>.jld2    (ComputeDissipation.jl)
```

The three scripts must be run in that order for a given point --
`ExcitedStateDMRG.jl` and `ComputeDissipation.jl` each load that point's
`Pss_Vdd<vDD_vT>_L<L>.h5` and will error out immediately (naming the file
they expected) if it isn't there yet. `Figures/AssembleTauAndDissipation.jl`
is the last step: once all 16 points are done, it reads every `mu1_*.jld2`
and `Dissipation_*.jld2`, derives $\langle\tau_{\rm err}\rangle = 2/\mu_1$,
and writes `Figures/"Figure 5 Data"/tau_err.txt` and `dissipation.txt` --
see that script's header comment for a caveat on where the $2/\mu_1$ formula
came from, and `Figures/"Figure 5 Data"/README.txt` for the exact commands.

## Per-figure regeneration chain

The path to creating the figure from scratch, data and all:

**Figure 1** (single-unit-scale steady-state density, marginalized over the
two right-most nodes of the L=3 reference chain):

```
Solvers/GetL3EigenvectorsWithGMRES/SparseSolve.jl   (runs the L=3 GMRES solve, marginalizes
                                              Pss over v2,v3, writes the raw export)
  -> Figures/Figure_1/Figure 1 Data/Figure1AxisTicks.txt, Figure1HeatMap.txt
  -> Figures/Figure_1/prep_fig1.py            (deterministic; no styling decisions)
  -> Figures/Figure_1/fig1_heatmap.dat, fig1_marginal.dat
  -> Figures/Figure_1/fig1.tex  ->  fig1.pdf
```

**Figure 2** (L=1 eigenbasis vectors, the five basis functions retained at the
production basis size):

```
Solvers/GetEigenvectorsWithDMRG/BasisTransform.jl (BuildU())  <-  called from ExcitedStateDMRG.jl
  -> Figures/Figure_2/Figure 2 Data/BasisVectorsXAxis.txt, BasisVectorsYAxis.txt
  -> Figures/Figure_2/prep_fig2.py
  -> Figures/Figure_2/fig2_basis.dat
  -> Figures/Figure_2/fig2.tex  ->  fig2.pdf
```


**Figure 3** (the paper's convergence argument: L2 error of the reduced-basis
DMRG eigenvectors, as a function of retained basis size $d$ and bond dimension
$\chi$, against the exact L=3 reference):

```
Solvers/GetL3EigenvectorsWithGMRES/L3_CMOS_System.jl + Build_Sparse_RateMatrix.jl + SparseSolve.jl
  -> Figures/Figure_3/Figure 3 Data/GMRES_Pss_L3.h5, GMRES_Prelax_L3.h5   (exact L=3 reference)
Figures/Figure_3/CompareGMRESwithTNsL3.jl   (sweeps d=2,4,...,64 and chi=10,...,90;
                                              see the design notes at the top of that
                                              file for exactly how "exact enumeration"
                                              and bond-dimension truncation are done)
  -> Figures/Figure_3/Figure 3 Data/L2err_Pss.txt, L2err_Prelax.txt
  -> Figures/Figure_3/prep_fig3.py
  -> Figures/Figure_3/fig3_pss.dat, fig3_prelax.dat, fig3_pss_c8.dat, fig3_prelax_c8.dat
  -> Figures/Figure_3/fig3.tex  ->  fig3.pdf
```

`CompareGMRESwithTNsL3.jl` needs `Figure 3 Data/GMRES_Pss_L3.h5` and
`GMRES_Prelax_L3.h5` to exist first, so run `SparseSolve.jl` before it.
Those two files are ~137 MB each and are intentionally not committed to
this repository -- they are reproducible raw data (`SparseSolve.jl`
regenerates them deterministically) rather than something that needs to
ship with the deposit; only the small `L2err_Pss.txt`/`L2err_Prelax.txt`
summary that Figure 3 actually plots is committed.

**Figures 4 and 5** (Arrhenius scaling and the main tau/dissipation tradeoff
result):

```
Solvers/GetEigenvectorsWithDMRG/SteadyStateDMRG.jl  <L> <vDD_vT>    (all 16 production points)
Solvers/GetEigenvectorsWithDMRG/ExcitedStateDMRG.jl <L> <vDD_vT>    (same point, after the above)
Solvers/ComputeDissipation/ComputeDissipation.jl    <L> <vDD_vT>    (same point, after the above)
  -> Solvers/GetEigenvectorsWithDMRG/Data/mu1_*.jld2, Solvers/ComputeDissipation/Data/Dissipation_*.jld2
Figures/AssembleTauAndDissipation.jl        (all 16 points must be present)
  -> Figures/"Figure 5 Data"/tau_err.txt, dissipation.txt
  -> Figures/Figure_4/prep_fig4.py  -> fig4_points.dat, fig4_fits.dat, fig4_r2.tex  -> fig4.tex -> fig4.pdf
  -> Figures/Figure_5/prep_fig5.py  -> fig5_tau.dat, fig5_dissipation.dat, fig5_tradeoff.dat, fig5_ref.dat -> fig5.tex -> fig5.pdf
```

## Software

Julia (this repo was developed against Julia 1.10) with `ITensors.jl`,
`ITensorMPS.jl`, `KrylovKit.jl`, `Arpack.jl`, `HDF5.jl`, `JLD2.jl`, and
`DelimitedFiles` (standard library). Figure prep scripts need Python 3 with
`numpy`. No version-pinned `Project.toml`/`Manifest.toml` is committed yet --
see "Known gaps."
