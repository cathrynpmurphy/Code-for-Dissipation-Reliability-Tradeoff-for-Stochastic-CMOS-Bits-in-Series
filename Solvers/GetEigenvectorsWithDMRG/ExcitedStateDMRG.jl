using ITensors, ITensorMPS, LinearAlgebra, JLD2, HDF5, DelimitedFiles

include("BuildMPO.jl") # Auxillary code for building the master equation as an MPO with ITensors AutoMPO()
include("MPOoperators.jl") # Auxillary code defining the operators
include("BasisTransform.jl") # Auxillary code for functions that perform the change of basis

#######################################
# Define the parameters of our circuit:
#######################################
#
# L and vDD_vT are read from the command line, exactly as in SteadyStateDMRG.jl (see the
# comment at the top of that file for the production shell-loop pattern) -- run this
# AFTER SteadyStateDMRG.jl has produced Data/Pss_Vdd<vDD_vT>_L<L>.h5 for the SAME (L,
# vDD_vT) point, since that's the file loaded below. Previously this file set its own
# vDD_vT = 1.1, independent of (and inconsistent with) SteadyStateDMRG.jl's vDD_vT = 1.2
# default -- i.e. the excited-state MPO used to be built at different physics than the
# steady-state penalty MPS loaded into it. Defaults are now shared (1.2) so a no-argument
# run of either file is self-consistent; please confirm 1.2 (not 1.1) is the value you
# actually want as the fallback.

# Parameters of each unit:
ve_vT =  0.1 #v_e / V_T -> the size of the discrete voltage step compared to the thermal voltage
nval = 1 #the slope factor n that parameterizes the transistors
vDD_vT = length(ARGS) >= 2 ? parse(Float64, ARGS[2]) : 1.2 #V_DD / V_T -> the size of the drain voltage compared to the thermal voltage
L = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 7 #The number of units in our daisy chain

###############################################
# Define the parameters of your Hilbert space for your occupation basis #
###############################################

# Define the limit of your finite state space (must match SteadyStateDMRG.jl's M, since
# that's what determines the shape of the Pss MPS being loaded below):
M = 32 # The largest number of (positive) discrete voltage units that your space spans
Mvals = [(M*(-1)):M;] # The voltage at each node v_i can take on the following discrete voltage values: [(M*(-1)):M;] .* ve_vT
n = length(Mvals) # Your original physical dimension

# Define the ITensors Hilbert space for the voltages at node i:
ITensors.space(::SiteType"Fock") = n
#Create the site indices for your MPO in your occupation basis with L units:
sitesOccBasis = siteinds("Fock",(L+1))

#########################################################################
## Import the Steady State Distribution computed by SteadyStateDMRG.jl ##
#########################################################################

# Import the steady state at this SAME (L, vDD_vT) point: assumes MPS form of the vector
# Pss(v0,v1...vL) in the occupation basis.
PathToSteadyState = "Data/"
PssFile = PathToSteadyState*"Pss_Vdd$(vDD_vT)_L$(L).h5"
isfile(PssFile) || error("$PssFile not found -- run SteadyStateDMRG.jl $L $vDD_vT first.")
Pss = h5open(PssFile, "r") do f
    read(f, "Pss", MPS)
end
# Defensive check: the loaded MPS should actually match the (L, M) this run just built --
# catches a stale/mismatched file slipping in under the expected name.
@assert length(siteinds(Pss)) == L+1 "Pss.h5's site count doesn't match L=$L"
@assert dim(siteinds(Pss)[1]) == n "Pss.h5's physical dimension doesn't match M=$M"

#################################################################
# Build your master equation as a Matrix Product Operator (MPO) #
#################################################################

H = MPO_Mechanisms(sitesOccBasis,"Series",L)
	# MPO_Mechanisms() -> function from "BuildMPO.jl"

##############################################################################
# Transform your MPO from the occupation basis to the L=1 eigenvector basis: #
##############################################################################

# Create the site indices for the L=1 system; used at every warm-start stage
# below to rebuild the linear map U for that stage's number of basis vectors.
sites1site = siteinds("Fock",2) #Create the site indices for the L=1 system
# We also also build it here to save out the first five vectors to create the data used in Figure 2
Ufig2 = BuildU(sites1site,M,5)
pathFig2 = "../../Figures/Figure_2/"

# Save Ufig2 out for Figures/Figure_2/prep_fig2.py, which reads
# "BasisVectorsXAxis.txt" (a length-61 axis) and "BasisVectorsYAxis.txt" (a
# 5x61 matrix, one row per retained basis vector f0..f4).
#
# Ufig2 itself depends on vDD_vT (BuildU's L=1 MPO uses the U+/U-/D+/D-
# operators, which do), so now that this file is driven across 16 different
# (L, vDD_vT) production points, only the canonical/default point actually
# writes Figure 2's export -- otherwise the last point run in the sweep would
# silently overwrite it with basis vectors from an arbitrary parameter choice.
if L == 7 && vDD_vT == 1.2
    mkpath(pathFig2*"Figure 2 Data")
    axisFull = Mvals .* ve_vT              # -3.2:0.1:3.2, 65 points
    WINDOW = 3:63                           # -3.0:0.1:3.0, 61 points (see note above)
    open(pathFig2*"Figure 2 Data/BasisVectorsXAxis.txt", "w") do io
        writedlm(io, axisFull[WINDOW])
    end
    open(pathFig2*"Figure 2 Data/BasisVectorsYAxis.txt", "w") do io
        writedlm(io, real.(Ufig2[:,WINDOW]))
    end
else
    println("Skipping Figure 2 export for (L=$L, vDD_vT=$vDD_vT) -- only written at the canonical L=7, vDD_vT=1.2 point.")
end

################################################################
## Create your initial state vector |P> as a matrix product state (MPS):
################################################################

#Create a uniform MPS in the occupation basis. This is the state that seeds
#the very first warm-start stage below.
Psi = MPS(sitesOccBasis)
nvec = ones(n) ./ n #uniform vector with length of physical dimension
for i=1:(L+1) #Loop over your nodes i
    Psi[i] = ITensor(nvec,sitesOccBasis[i]) #Initialize your state vector |P> as a uniform product state
end

##########################################################################
## WARM-START EXCITED STATE DMRG (Sec. IV D of the paper) ################
##########################################################################
# As in SteadyStateDMRG.jl, DMRG is converged in stages rather than run
# once at full resolution: we begin from the uniform MPS above with a
# small maximal bond dimension (chi = MD) and few retained single-unit
# basis vectors (d = nbasis), then increase both together as the MPS
# converges, ending at the production values (d = 32, chi <= 100) used to
# report results in Sec. V. Because the steady state |Pss> is used at
# every stage as the DMRG penalty state, it must be re-transformed
# (Pss2 = TransformPsi(Pss,...)) into each stage's basis before that
# stage's DMRG calls, using that stage's U -- reusing a Pss2 built in a
# different stage's basis would silently penalize the wrong state.
#
# Each row of the schedule is one stage: [nbasis (d), maxdim (chi), ncycles].
# As in SteadyStateDMRG.jl, we step d up by 2 basis vectors per stage,
# from the starting point (d = 12) to the production endpoint (d = 32),
# giving 11 stages. Each stage runs ~540 single-site sweeps (matching the
# paper's total of ~5400 sweeps for the hardest case, L = 7,
# Vdd/VT = 1.4, spread evenly over the 11 stages); with StepNum DMRG
# iterations per cycle below, that's ncyclesPerStage = 540 / StepNum
# cycles per stage. The bond dimension chi is ramped in step with d, from
# a small starting value up to the production cap (chi <= 100), per
# "both cutoffs are increased together" in Sec. IV D. Ramped 80 -> 100,
# identically to SteadyStateDMRG.jld2
StepNum=5 #Number of DMRG iterations per cycle
sweepsPerStage = 540 #Total single-site DMRG sweeps to run at each stage
ncyclesPerStage = Int(round(sweepsPerStage/StepNum)) #Number of DMRG-function calls ("cycles") per stage

nbasisSchedule = collect(12:2:32) #d: 12,14,16,...,32 -> 11 warm-start stages
nstages = length(nbasisSchedule)
maxdimSchedule = round.(Int,range(80,100,length=nstages)) #chi: ramped in step with d, from 80 up to the production cap of 100

warmStartSchedule = hcat(nbasisSchedule, maxdimSchedule, fill(ncyclesPerStage,nstages))

## DMRG Parameters shared by every stage: #####
#########################
cutoff=1e-20 #Truncation cutoff for singular values
olevel=1 #Readout level while DMRG runs
#########################

energies = Complex[] # A running record of the inner product value <Psi|H|Psi>/<Psi|Psi>, across every stage/cycle

# Created up front, not just before the final Prelax.h5 save below -- the per-cycle
# Energies_Prelax_*.jld2 checkpoint inside the loop below writes here starting on the
# very first cycle.
mkpath("Data")

for stage = 1:size(warmStartSchedule,1)

	nbasis = warmStartSchedule[stage,1] # The number of basis vectors at this stage (your new physical dimension)
	MD     = warmStartSchedule[stage,2] # The maximum bond dimension DMRG will permit before truncating at this stage
	ncycles= warmStartSchedule[stage,3] # Number of times we'll call the DMRG function at this stage

	println("Warm-start stage $stage: d = $nbasis basis vectors, max bond dim = $MD, $ncycles cycles.")

	# Create a linear map from the occupation basis to this stage's new basis
	U = BuildU(sites1site,M,nbasis) #Solve the L=1 system and turn its eigenvectors into a matrix
		#BuildU() -> Function from "BasisTransform.jl"

	# Define the ITensors Hilbert space for the new basis vectors at node i:
	ITensors.space(::SiteType"EigenBasis") = nbasis
	# Create the site indices for your MPO in your new basis with L units:
	sitesNewBasis = siteinds("EigenBasis",(L+1))

	# Transform the MPO, the best current guess for Psi, and the steady
	# state penalty vector, into this stage's basis:
	H2 = TransformH(H,sitesNewBasis,U)
		#TransformH() -> Function from "BasisTransform.jl"
	Psi2 = TransformPsi(Psi,sitesNewBasis,U)
		#TransformPsi() -> Function from "BasisTransform.jl"
	Pss2 = TransformPsi(Pss,sitesNewBasis,U)
		#TransformPsi() -> Function from "BasisTransform.jl"

	for t=1:ncycles
		#Call the ITensor DMRG Function with a penalty against [Pss2]
		Henergy, PsiOut = dmrg(((-1)*H2),[Pss2],Psi2; nsweeps=StepNum,cutoff=cutoff,maxdim=MD,outputlevel=olevel,ishermitian=false)
		#Save the outputs:
		push!(energies,Henergy) #This should converge to first nonzero eigenvalue
		h5open("PsiOut.h5", "w") do f
	    	write(f, "PsiOut", PsiOut)
		end
		#Update your state vector MPS:
		Psi2 .= PsiOut
	end

	#Bring this stage's converged MPS back to the occupation basis so it can
	#seed the next (higher-resolution) stage as a warm start:
	Psi = TransformPsiBack(Psi2,sitesOccBasis,U)
		#TransformPsiBack() -> Function from "BasisTransform.jl"

end

#Psi is now the MPS from the final (production-basis) warm-start stage,
#already transformed back into the occupation basis, i.e. the first excited
#state eigenvector Prelax(v0,v1,...vL)
Prelax = Psi
h5open("Data/Prelax_Vdd$(vDD_vT)_L$(L).h5", "w") do f
    	write(f, "Prelax", Prelax)
end

### The eigenvalue associated with |P_relax>, the slowest switching timescale is the converged excited state DMRG energy:
mu1 = energies[length(energies)]
# Named "mu1" (ascii), not "μ1" as before -- unicode filenames are an avoidable
# cross-platform/tooling risk for a Zenodo deposit.
save_object("Data/mu1_Vdd$(vDD_vT)_L$(L).jld2",mu1)
println("The spectral gap of the circuit master equation is $mu1.")
