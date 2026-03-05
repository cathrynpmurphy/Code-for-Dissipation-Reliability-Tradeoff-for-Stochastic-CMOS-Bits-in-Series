Main Code Files:
- SteadyStateDMRG.jl -> Computes the steady state eigenvector for a CMOS circuit with parameters {k} and L repeated units
- ExcitedStateDMRG.jl -> Uses the penalty method to penalize the steady state and compute the first excited state eigenvector, corresponding to the slowest relaxation mode

Auxillary Code Files:
- MPOoperators.jl -> Defines the second quantized operators that will be utilized to build the matrix product operators
- BuildMPO.jl -> Constructs the matrix product operator (MPO) representation of the CMOS master equation for a particular parameter set and number of CMOS units
- BasisTransform.jl -> A set of functions which are used to transform the dynamics (MPO) and state vector (MPS) between the Fock basis and the L=1 eigenvector basis