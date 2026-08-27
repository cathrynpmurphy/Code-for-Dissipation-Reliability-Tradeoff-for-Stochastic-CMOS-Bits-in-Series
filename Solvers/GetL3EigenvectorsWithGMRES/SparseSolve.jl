using KrylovKit
using ITensors
using JLD2
using LinearAlgebra
using Arpack
using HDF5
using DelimitedFiles

flush(stdout)

## Assemble the sparse CMOS rate matrix operator using the auxillary code:
include("Build_Sparse_RateMatrix.jl")
include("L3_CMOS_System.jl")
    # Defines H_tot, the L=3 sparse rate matrix, along with n, M, Mvals, ve_vT
    # (needed below to write out the Figure 1 axis values).

## Call the GMRES solver to get the two eigenvectors associated with the two eigenvalues closest to zero:
vals,vecs  = KrylovKit.eigsolve(H_tot, 2, :LR; krylovdim=50, maxiter=100, verbosity=2)

Pss = reshape(vecs[1],(n,n,n,n)) # Zero eigenvector
Prelax = reshape(vecs[2],(n,n,n,n)) # Eigenvector with smallest nonzero eigenvalues

## Save out the full distributions to compare with the distributions truncated using Tensor networks and the basis transformations:
# Save the data both to the Figure 3 path and to this folder, so the deposit's Figure 3
# chain (Figures/Figure_3/CompareGMRESwithTNsL3.jl -> prep_fig3.py -> fig3.tex) is
# self-contained and doesn't reach outside this repository.
pathFig3 = "../../Figures/Figure_3/"
mkpath(pathFig3*"Figure 3 Data")

h5open(pathFig3*"Figure 3 Data/GMRES_Pss_L3.h5", "w") do f
    write(f, "Pss", Pss)
end
h5open(pathFig3*"Figure 3 Data/GMRES_Prelax_L3.h5", "w") do f
    write(f, "Prelax", Prelax)
end

## Marginalize Pss in terms of v0 and v1; this is used to generate the data in Figure 1
Pss_v0v1 = dropdims(sum(Pss, dims=(3,4)), dims=(3,4))
Pss_v0v1 /= sum(Pss_v0v1) # Normalize so that the marginalized distribution sums to ones
axisVals = [(-M*ve_vT):ve_vT:(M*ve_vT);]

## Save out the marginalized distribution and axis values as data for the heat map in Figure 1
# Save the data both to the Figure 1 path and to this folder. Written as plain
# tab-separated text (via DelimitedFiles.writedlm) so Figures/Figure_1/prep_fig1.py's
# np.loadtxt(...) reads them directly, matching the format PROMPT_fig1.md documents:
# Figure1AxisTicks.txt (65 values) and Figure1HeatMap.txt (65x65 matrix).
pathFig1 = "../../Figures/Figure_1/"
mkpath(pathFig1*"Figure 1 Data")

writedlm(pathFig1*"Figure 1 Data/Figure1AxisTicks.txt", axisVals, '\t')
# real(...) because Pss_v0v1 is complex (KrylovKit returns complex eigenvectors even
# though the physical steady-state distribution is real to numerical precision); the
# imaginary part should be ~0 and prep_fig1.py's own assertions (normalisation, no
# large negative entries) are the check that this projection was legitimate.
writedlm(pathFig1*"Figure 1 Data/Figure1HeatMap.txt", real.(Pss_v0v1), '\t')

println("Wrote Pss.h5, Prelax.h5, and the Figure 1/Figure 3 raw data files.")
