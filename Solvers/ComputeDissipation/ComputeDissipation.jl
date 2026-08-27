using ITensors
using LinearAlgebra
using Pkg
using HDF5
using JLD2

include("BuildDissipationMPO.jl")
include("DissipationOperators.jl")
include("DissipationFunctions.jl")

#######################################
# Define the parameters of our circuit:
#######################################
#
# L and vDD_vT are read from the command line, exactly as in SteadyStateDMRG.jl /
# ExcitedStateDMRG.jl -- run this AFTER SteadyStateDMRG.jl has produced
# ../GetEigenvectorsWithDMRG/Data/Pss_Vdd<vDD_vT>_L<L>.h5 for the SAME (L, vDD_vT) point.
# Previously this file set its own vDD_vT = 1.1, independent of (and, before this pass,
# inconsistent with) SteadyStateDMRG.jl's default -- the dissipation MPO would have been
# built at different physics than the Pss it was applied to. Defaults are now shared
# (1.2) with the other two files.

# Parameters of each unit:
ve_vT =  0.1 #v_e / V_T -> the size of the discrete voltage step compared to the thermal voltage
nval = 1 #the slope factor n that parameterizes the transistors
vDD_vT = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.2 #V_DD / V_T -> the size of the drain voltage compared to the thermal voltage
L = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 7 #The number of units in our daisy chain

###############################################
# Define the parameters of your Hilbert space for your occupation basis #
###############################################

# Must match SteadyStateDMRG.jl's M, since that's what determines the shape of the Pss
# MPS being loaded below:
M = 32 # The largest number of (positive) discrete voltage units that your space spans
Mvals = [(M*(-1)):M;] # The voltage at each node v_i can take on the following discrete voltage values: [(M*(-1)):M;] .* ve_vT
n = length(Mvals) # Your original physical dimension

# Define the ITensors Hilbert space for the voltages at node i:
ITensors.space(::SiteType"Fock") = n
sitesOccBasis = siteinds("Fock",(L+1))

#########################################################################
## Import the Steady State Distribution computed by SteadyStateDMRG.jl ##
#########################################################################

# Import the steady state at this SAME (L, vDD_vT) point: assumes MPS form of the vector
# Pss(v0,v1...vL) in the occupation basis.
PathToSteadyState = "../GetEigenvectorsWithDMRG/Data/"
PssFile = PathToSteadyState*"Pss_Vdd$(vDD_vT)_L$(L).h5"
isfile(PssFile) || error("$PssFile not found -- run SteadyStateDMRG.jl $L $vDD_vT first.")
Pss = h5open(PssFile, "r") do f
    read(f, "Pss", MPS)
end
@assert length(siteinds(Pss)) == L+1 "Pss.h5's site count doesn't match L=$L"
@assert dim(siteinds(Pss)[1]) == n "Pss.h5's physical dimension doesn't match M=$M"

# Normalize the steady state distribution.
# NOTE: this used to be `Pss = NormalizePss(Pss)`, which assigns NormalizePss's whole
# (PssNorm, Nval) tuple return to Pss -- silently turning Pss into a 2-tuple instead of
# an MPS, which would have broken `siteinds(Pss)` right below it. Fixed to destructure
# both return values.
Pss, PssNorm = NormalizePss(Pss)
  # NormalizePss() -> function from "DissipationFunctions.jl"

#########################################
# Create the MPO to compute dissipation #
#########################################

# This operator (built as an MPO) is constructed such that <Q̇> = <1|Hdiss|Pss>
Hdiss = Dissipation_MPO(sitesOccBasis,"Series",L)

# Loop over the tensor network sites to compute the inner product <1|Hdiss|Pss>:
# Pss[i] -> factorized ket vector |Pss(vi)>
# Hdiss[i] -> factorized operator Hdiss(vi,vi')
# ITensor(ones(1,n),sitesOccBasis[i]) -> Bra vector <1| to sum over all values of vi
#
# NOTE: this loop runs over sites 1:(L+1) -- there are L+1 sites for an L-unit chain
# (v0..vL). The previous version of this file derived its own `L` as
# `length(sitesOccBasis)`, i.e. the SITE COUNT (L_chain+1), and used that same value both
# here (where site-count is in fact what's needed, so this loop was correct) AND as the
# chain-length argument to Dissipation_MPO two lines up (where it is NOT what's needed --
# Dissipation_MPO's "Series" branch loops site indices up to L+1, so passing L_chain+1
# there asks it to reach one site past the end of sitesOccBasis and would have errored).
# Now that `L` consistently means chain length everywhere in this file (matching
# SteadyStateDMRG.jl/ExcitedStateDMRG.jl), this loop's bound is L+1, not L.
Qtot = Pss[1]*Hdiss[1]*ITensor(ones(1,n),sitesOccBasis[1])
for i=2:(L+1)
   global Qtot
   Qtot *= Pss[i]*Hdiss[i]*ITensor(ones(1,n),sitesOccBasis[i])
end

Qtot = scalar(Qtot)
mkpath("Data")
save_object("Data/Dissipation_Vdd$(vDD_vT)_L$(L).jld2",Qtot)

println("The Dissipation Rate at (L=$L, vDD_vT=$vDD_vT) is $Qtot.")
