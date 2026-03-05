This Folder Contains: 

- GetEigenvectorsWithDMRG -> A folder containing the necessary code to obtain the steady state eigenvector and first nonzero eigenvector by:
    - Constructing the CMOS circuit dynamics as a matrix product operator (MPO)
    - Applying a basis transformation to the dynamics and the probability state space
    - Using DMRG to compute the steady state eigenvector
    - Using DMRG with a penalty method to compute the first non-zero eigenvector, which corresponds to the slowest dynamical mode and the spectral gap

- ComputeDissipation -> A folder containing the code necessary to compute the dissipation associated with a circuit given its steady state eigenvector by:
    - Constructing an operator as an MPO whose steady state expectation value is the dissipation rate
    - Applying the dissipation MPO sequentially to the steady state MPS and computing the inner product