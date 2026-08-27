Main Code Files:
- SparseSolve.jl -> Computes the steady state and first excited eigenvector with GMRES at L=3
  (the exact reference Figures/Figure_3/CompareGMRESwithTNsL3.jl checks the tensor-network
  result against). Also writes the Figure 1 and Figure 3 raw data exports -- see the
  top-level README.md for the full per-figure chain.

Auxillary Code Files:
- L3_CMOS_System.jl -> Encodes the parameter information about your CMOS chain (L=3) and
  builds the per-site operators (x+, x-, U+, U-, D+, D-, upperProj, lowerProj), named to
  match Solvers/GetEigenvectorsWithDMRG/MPOoperators.jl
- Build_Sparse_RateMatrix.jl -> Assembles those operators into the sparse rate matrix for
  a chain of L units in series, mirroring Solvers/GetEigenvectorsWithDMRG/BuildMPO.jl term-by-term