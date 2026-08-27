using SparseArrays

#######################################
# Define the parameters of our circuit:
#######################################

# Parameters of each unit -- these must match Solvers/GetEigenvectorsWithDMRG/SteadyStateDMRG.jl
# and ExcitedStateDMRG.jl for the L=3 comparison in Figure 3 to be a fair test of the
# tensor-network method against this exact/sparse reference:
ve_vT =  0.1 #v_e / V_T -> the size of the discrete voltage step compared to the thermal voltage
nval = 1 #the slope factor n that parameterizes the transistors
vDD_vT = 1.2 #V_DD / V_T -> the size of the drain voltage compared to the thermal voltage

# The number of units in our daisy chain:
# NOTE: this is the L=3 reference solve that Figure 3 validates the reduced-basis DMRG
# result against (see Figures/Figure_3/PROMPT_fig345.md: "Two stacked panels at L=3").
# L=3 keeps the exact state space small enough (n^4, below) for an iterative sparse
# eigensolver; it is NOT the L=7 production chain length used to report the paper's main
# results in Solvers/GetEigenvectorsWithDMRG. 
L=3 #The number of units in our chain (matches the "L3" in this folder's name)

###############################################
# Define the parameters of your Hilbert space for your occupation basis #
###############################################

# Define the limit of your finite state space:
M = 32 # The largest number of (positive) discrete voltage units that your space spans
Mvals = [(M*(-1)):M;] # The voltage at each node v_i can take on the following discrete voltage values: [(M*(-1)):M;] .* ve_vT
n = length(Mvals) # Your original physical dimension

#####################################################################################
# Build the sparse operators, using the SAME operator names as the tensor-network   #
# code (Solvers/GetEigenvectorsWithDMRG/MPOoperators.jl) so the two independent solvers can #
# be compared term-by-term. These are the matrix (not ITensor) versions of the same #
# five physical operators: x+/x- (raising/lowering) and U+/U-/D+/D- (the PMOS/NMOS  #
# uphill/downhill diagonal rates).                                                  #
#####################################################################################

## Raising / lowering operators ("x+" / "x-" in MPOoperators.jl) ##
Xminus = spzeros(n,n) # lowering operator, matches ::OpName"x-"
for a=1:(2*M)
    Xminus[a,a+1] = 1
end

Xplus = spzeros(n,n) # raising operator, matches ::OpName"x+"
for a=1:(2*M)
    Xplus[a+1,a] = 1
end

## Operator that projects out the upper state to ensure probability conservation ##
upperProj = spzeros(n,n)
for a=1:(2*M)
    upperProj[a,a] = 1
end

## Operator that projects out the lowest state to ensure probability conservation ##
lowerProj = spzeros(n,n)
for a=1:(2*M)
    lowerProj[a+1,a+1] = 1
end

## Diagonal rate operators, identical formulas to MPOoperators.jl's U+/U-/D+/D- ##
Uplus = spzeros(n,n) # "U+": Uphill associated with PMOS
Dplus = spzeros(n,n) # "D+": Downhill associated with PMOS
Uminus = spzeros(n,n) # "U-": Uphill associated with NMOS
Dminus = spzeros(n,n) # "D-": Downhill associated with NMOS
for a=1:(2*M+1)
    mval = Mvals[a]
    Uplus[a,a] = exp(vDD_vT/nval)*exp(-1*mval*ve_vT/nval)
    Dplus[a,a] = exp((-1/2)*ve_vT-vDD_vT)*exp(mval*ve_vT)
    Uminus[a,a] = exp(vDD_vT/nval)*exp(mval*ve_vT/nval)
    Dminus[a,a] = exp((-1/2)*ve_vT-vDD_vT)*exp(-1*mval*ve_vT)
end

# Identity operator
ID = spzeros(n,n)
for a=1:(2*M + 1)
    ID[a,a] = 1
end

#-------------------------------------------------------------------------------
H_tot = H_Mechanisms(L,n,Xplus,Xminus,Uplus,Uminus,Dplus,Dminus,upperProj,lowerProj,ID,"Series") # Call the function to create the sparse rate matrix operator
#-------------------------------------------------------------------------------
