using ITensors, ITensorMPS, LinearAlgebra, HDF5, JLD2

include("BuildMPO.jl") # Auxillary code for building the master equation as an MPO with ITensors AutoMPO()
include("MPOoperators.jl") # Auxillary code defining the operators
include("BasisTransform.jl") # Auxillary code for functions that perform the change of basis

#######################################
# Define the parameters of our circuit:
#######################################

# Parameters of each unit:
ve_vT =  0.1 #v_e / V_T -> the size of the discrete voltage step compared to the thermal voltage
nval = 1 #the slope factor n that parameterizes the transistors
vDD_vT = 1.1 #V_DD / V_T -> the size of the drain voltage compared to the thermal voltage

# The number of units in our daisy chain:
L=7 #The number of units in our chain

###############################################
# Define the parameters of your Hilbert space for your occupation basis #
###############################################

# Define the limit of your finite state space:
M = 32 # The largest number of (positive) discrete voltage units that your space spans
Mvals = [(M*(-1)):M;] # The voltage at each node v_i can take on the following discrete voltage values: [(M*(-1)):M;] .* ve_vT
n = length(Mvals) # Your original physical dimension

# Define the ITensors Hilbert space for the voltages at node i:
ITensors.space(::SiteType"Fock") = n
#Create the site indices for your MPO in your occupation basis with L units:
sitesOccBasis = siteinds("Fock",(L+1))

#################################################################
# Build your master equation as a Matrix Product Operator (MPO) #
#################################################################

H = MPO_Mechanisms(sitesOccBasis,"Series",L)
	# MPO_Mechanisms() -> function from "BuildMPO.jl"

##############################################################################
# Transform your MPO from the occupation basis to the L=1 eigenvector basis: #
##############################################################################

# Decide how many basis vectors in your transformed basis you would like to retain
nbasis = 12 # The number of basis vectors (your new physical dimension)

# Create a linear map from the occupation basis to the new basis
sites1site = siteinds("Fock",2) #Create the site indices for the L=1 system
U = BuildU(sites1site,M,nbasis) #Solve the L=1 system and turn its eigenvectors into a matrix
	#BuildU() -> Function from "BasisTransform.jl"

# Define the ITensors Hilbert space for the new basis vectors at node i:
ITensors.space(::SiteType"EigenBasis") = nbasis
# Create the site indices for your MPO in your new basis with L units:
sitesNewBasis = siteinds("EigenBasis",(L+1))

# Transform the old MPO into the MPO in the new basis:
H2 = TransformH(H,sitesNewBasis,U)
	#TransformH() -> Function from "BasisTransform.jl"

################################################################
## Create your state vector |P> as a matrix product state (MPS):
################################################################

#Create a uniform MPS
Psi = MPS(sitesOccBasis)
nvec = ones(n) ./ n #uniform vector with length of physical dimension
for i=1:(L+1) #Loop over your nodes i
    Psi[i] = ITensor(nvec,sitesOccBasis[i]) #Initialize your state vector |P> as a product state of the L=1 steady state eigenvector
end

# Transform the old MPO into the MPS in the new basis:
Psi2 = TransformPsi(Psi,sitesNewBasis,U)
	#TransformPsi() -> Function from "BasisTransform.jl"

############################
## BEGINNING DMRG ##########
############################

## DMRG Parameters: #####
#########################
ncycles=200 # Number of times we'll call the DMRG function
cutoff=1e-20 #Truncation cutoff for singular values
MD=100 #The maximum bond dimension DMRG will permit before truncating
olevel=1 #Readout level while DMRG runs
StepNum=5 #Number of DMRG iterations per cycle
#########################

energies = zeros(Complex,ncycles) # A running record of the inner product value <Psi|H|Psi>/<Psi|Psi> as DMRG iterates on Psi

for t=1:ncycles 
	#Call the ITensor DMRG Function:
	Henergy, PsiOut = dmrg(((-1)*H2),Psi2; nsweeps=StepNum,cutoff=cutoff,maxdim=MD,outputlevel=olevel,ishermitian=false)	
	#Save the outputs:
	energies[t] = Henergy # This should converge to zero
	save_object("Energies.jld2",energies[1:t])
	h5open("PsiOut.h5", "w") do f
    	write(f, "PsiOut", PsiOut)
	end
	#Update your state vector MPS:
	Psi2 .= PsiOut 
end

#Take the MPS converged with DMRG and transform it back to the original occupation basis
# in order to extract Pss(v0,v1,...vL)
Pss = TransformPsiBack(Psi2,sitesOccBasis,U)
# Save the steady state MPS as a .h5 file
h5open("Pss.h5", "w") do f
    write(f, "Pss", Pss)
end
