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

# Parameters of each unit:
ve_vT =  0.1 #v_e / V_T -> the size of the discrete voltage step compared to the thermal voltage
nval = 1 #the slope factor n that parameterizes the transistors
vDD_vT = 1.1 #V_DD / V_T -> the size of the drain voltage compared to the thermal voltage

#########################################################################
## Import the Steady State Distribution (should come from normal DMRG) ##
#########################################################################

# Import the steady state: assumes MPS form of the vector Pss(v0,v1...vL) in the occupation basis
PathToSteadyState = "/Path/To/SteadyState/"
Pss = h5open("$PathToSteadyState"*"/Pss.h5", "r") do f
    read(f, "Pss", MPS)
end

# Normalize the steady state distribution
Pss = NormalizePss(Pss)
  # NormalizePss() -> function from "Functions.jl"

#####################################################################################
### Learn the parameters of your ITensors Hilbert Space from the Steady State MPS: ##
#####################################################################################

# Read the ITensors site indices from your steady state MPS
sitesOccBasis = siteinds(Pss)
n = dim(sitesOccBasis[1]) #Your physical dimension in the occupation basis, which we get from the steady state MPS
M = floor(Int,((n-1)/2)) # The largest number of (positive) discrete voltage units that your space spans
Mvals = [(M*(-1)):M;] # The voltage at each node v_i can take on the following discrete voltage values: [(M*(-1)):M;] .* ve_vT

# The number of units in our daisy chain:
L = length(sitesOccBasis) #The number of units in our chain

# Define the ITensors Hilbert space for the voltages at node i:
ITensors.space(::SiteType"Fock") = n

#########################################
# Create the MPO to compute dissipation #
#########################################

# This operator (built as an MPO) is constructed such that <Q̇> = <1|Hdiss|Pss>
Hdiss = Dissipation_MPO(sitesOccBasis,"Series",L)

# Loop over the tensor network sites to compute the inner product <1|Hdiss|Pss>:
# Pss[i] -> factorized ket vector |Pss(vi)>
# Hdiss[i] -> factorized operator Hdiss(vi,vi')
# ITensor(ones(1,n),sitesOccBasis[i]) -> Bra vector <1| to sum over all values of vi
Qtot = Pss[1]*Hdiss[1]*ITensor(ones(1,n),sitesOccBasis[1])
for i=2:L 
   global Qtot
   Qtot *= Pss[i]*Hdiss[i]*ITensor(ones(1,n),sitesOccBasis[i])
end

Qtot = scalar(Qtot)
save_object("Dissipation.jld2",Qtot)

println("The Dissipation Rate is $Qtot.")
