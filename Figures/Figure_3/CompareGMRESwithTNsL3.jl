using ITensors, ITensorMPS
using HDF5
using DelimitedFiles
using LinearAlgebra

# Get the supplementary code that lets you perform a basis transformation, and the code
# that builds the L=3 master equation as an MPO (local copies, matching the pattern
# already used for BasisTransform.jl in this folder:
include("BuildMPO.jl")       # MPO_Mechanisms()
include("MPOoperators.jl")   # the "Fock"-site x+/x-/U+/U-/D+/D-/upperProj/lowerProj ops
include("BasisTransform.jl") # BuildU(), TransformH(), TransformPsi(), TransformPsiBack()

NSWEEPS_EXACT = 20  # generous sweep count; see note 1 above -- with cutoff=0 and
                     # maxdim=d^2 there is no truncation to converge past, only the
                     # variational optimization itself, which is fast for 4 sites.

#######################################
# Define the parameters of our circuit -- MUST match Solvers/GetL3EigenvectorsWithGMRES/L3_CMOS_System.jl
#######################################
ve_vT = 0.1
nval = 1
vDD_vT = 1.2
L = 3

M = 32
Mvals = [(M*(-1)):M;]
n = length(Mvals) # n = 65

ITensors.space(::SiteType"Fock") = n
sitesOccBasis = siteinds("Fock",(L+1))

# Build the L=3 master equation MPO in the occupation basis:
H = MPO_Mechanisms(sitesOccBasis,"Series",L)

# Site indices for the L=1 system, used by BuildU at every basis size below:
sites1site = siteinds("Fock",2)

################################################################
# Load the GMRES reference (Pss, Prelax), written by
# Solvers/GetL3EigenvectorsWithGMRES/SparseSolve.jl into this figure's local data folder.
################################################################
GMRESpath = "Figure 3 Data/"
Pss_gmres = h5open(GMRESpath*"GMRES_Pss_L3.h5", "r") do f
    read(f, "Pss")
end
Prelax_gmres = h5open(GMRESpath*"GMRES_Prelax_L3.h5", "r") do f
    read(f, "Prelax")
end

# L1-normalize both reference distributions (matches NormalizePss's convention: <1|P>=1):
Pss_gmres = real.(Pss_gmres) ./ sum(real.(Pss_gmres))
# Prelax is not a probability distribution (it integrates to ~0, it's the first excited
# eigenvector), so normalize it in the L2 sense instead, which is the natural gauge to
# compare eigenvectors that are only defined up to overall scale:
Prelax_gmres = real.(Prelax_gmres) ./ norm(real.(Prelax_gmres))

function mps_to_occbasis_array(psi, sitesOccBasis)
    # Contract an (L+1)-site MPS in the occupation basis down to a dense n x n x n x n
    # array, matching the shape GMRES's Pss/Prelax were reshaped into in SparseSolve.jl.
    full = psi[1]
    for i=2:length(psi)
        full *= psi[i]
    end
    return Array(full, sitesOccBasis...)
end

function l2_error_pss(psi, sitesOccBasis, ref)
    arr = mps_to_occbasis_array(psi, sitesOccBasis)
    arr = real.(arr) ./ sum(real.(arr))
    # sign/phase ambiguity doesn't apply to a normalized probability distribution
    # (it's non-negative by construction up to numerical noise), so no alignment needed
    return norm(arr .- ref)
end

function l2_error_prelax(psi, sitesOccBasis, ref)
    arr = mps_to_occbasis_array(psi, sitesOccBasis)
    arr = real.(arr) ./ norm(real.(arr))
    # align sign: flip if that reduces the distance to the reference
    if norm(arr .- ref) > norm(-arr .- ref)
        arr = -arr
    end
    return norm(arr .- ref)
end

########################################
# Sweep over basis size d = 2,4,...,64 #
########################################
basisSizes = collect(2:2:64)     # 32 values -- x-axis of Figure 3
bondDims   = collect(10:10:90)   # 9 values  -- y-axis of Figure 3

L2err_Pss    = zeros(ComplexF64, length(bondDims), length(basisSizes))
L2err_Prelax = zeros(ComplexF64, length(bondDims), length(basisSizes))

cutoff = 0.0
olevel = 1

# Seed state for every dmrg() call below: a uniform product state in the occupation
# basis, transformed into each basis size's new basis via TransformPsi -- the same
# seeding convention SteadyStateDMRG.jl/ExcitedStateDMRG.jl use for their warm start,
# reused here rather than inventing a separate seed convention for "EigenBasis" sites
# (which, unlike "Fock", has no named single-site states defined for MPS(sites,states)
# to look up).
PsiOccSeed = MPS(sitesOccBasis)
nvecSeed = ones(n) ./ n
for i=1:(L+1)
    PsiOccSeed[i] = ITensor(nvecSeed, sitesOccBasis[i])
end

for (jd, d) in enumerate(basisSizes)

    println("Basis size d = $d ($jd / $(length(basisSizes)))")

    # Build the transformation into this basis size's eigenbasis:
    U = BuildU(sites1site, M, d)

    ITensors.space(::SiteType"EigenBasis") = d
    sitesNewBasis = siteinds("EigenBasis",(L+1))

    H2 = TransformH(H, sitesNewBasis, U)

    maxdimExact = d^2  # unrestricted bond dimension for this basis size (see note 1)

    # ---- Pss: dominant (zero) eigenvector, plain dmrg minimizing -H2 -----------------
    Psi0 = TransformPsi(PsiOccSeed, sitesNewBasis, U)
    _, PssExact = dmrg(((-1)*H2), Psi0; nsweeps=NSWEEPS_EXACT, cutoff=cutoff,
                        maxdim=maxdimExact, outputlevel=olevel, ishermitian=false)

    # ---- Prelax: first-excited eigenvector, penalty dmrg against PssExact ------------
    Psi1 = TransformPsi(PsiOccSeed, sitesNewBasis, U)
    _, PrelaxExact = dmrg(((-1)*H2), [PssExact], Psi1; nsweeps=NSWEEPS_EXACT, cutoff=cutoff,
                           maxdim=maxdimExact, outputlevel=olevel, ishermitian=false)

    for (jchi, chi) in enumerate(bondDims)

        PssChi = deepcopy(PssExact)
        truncate!(PssChi; maxdim=chi, cutoff=0.0)
        PssChiOcc = TransformPsiBack(PssChi, sitesOccBasis, U)
        L2err_Pss[jchi,jd] = l2_error_pss(PssChiOcc, sitesOccBasis, Pss_gmres)

        PrelaxChi = deepcopy(PrelaxExact)
        truncate!(PrelaxChi; maxdim=chi, cutoff=0.0)
        PrelaxChiOcc = TransformPsiBack(PrelaxChi, sitesOccBasis, U)
        L2err_Prelax[jchi,jd] = l2_error_prelax(PrelaxChiOcc, sitesOccBasis, Prelax_gmres)

    end
end

########################################
# Save the raw L2 error grids.         #
########################################
# Written in the same "9x32 tab-separated Julia complex" format prep_fig3.py already
# expects (see its docstring): writedlm's default element formatting for Complex
# numbers matches Julia's string(::Complex) form ("a + bim"), which is what
# prep_fig3.py's load_complex() parses by splitting each cell on '+'.
mkpath("Figure 3 Data")
writedlm("Figure 3 Data/L2err_Pss.txt", L2err_Pss, '\t')
writedlm("Figure 3 Data/L2err_Prelax.txt", L2err_Prelax, '\t')

println("Wrote Figure 3 Data/L2err_Pss.txt and L2err_Prelax.txt ($(length(bondDims))x$(length(basisSizes))).")
println("Run prep_fig3.py next to turn these into fig3_pss.dat / fig3_prelax.dat and the 1e-8 contours.")
