###############################################
## Create the Linear Map from the L=1 System ##
###############################################

function BuildU(sites1site,M,nbasis)

    # Arguments: 
        # sites1site -> ITensor site indices for L=1 system
        # M -> Even Cap on Discrete Voltage Space
        # nbasis -> Number of basis vectors we would like to reserve
    # Returns:
	    # U -> A matrix whose rows are our new basis vectors as function of the discrete voltage values

    #Create the MPO for the L=1 Problem
    H1site = MPO_Mechanisms(sites1site,"SingleUnit",1)

    #Turn the MPO into a matrix
    Htens = H1site[1]*H1site[2]
    comb = combiner(sites1site)
    Htens *= comb
    Htens *= prime(comb,plev=0,1)
    Hmatrix = Matrix(Htens,inds(Htens)[2],inds(Htens)[1])

    #Solve for the eigenvectors of the L=1 MPO, which have the form P(v0,v1)
    eigset = eigen(Hmatrix)
    eigenvecs = eigset.vectors
    eigenvals = eigset.values

    #Define the original physical dimension
    n = 2*M + 1

    #Create a matrix whose columns will be our basis vectors
    U2 = zeros(Complex,n,nbasis)

    # Marginalize over one dimension to turn P(v0,v1) -> P(v0)
    for i=1:nbasis #Loop over all the eigenvectors
        vec = eigenvecs[:,(length(eigenvals)+1-i),:] #Extract the appropriate eigenvector
        mat = reshape(vec,(n,n)) #Reshape the eigenvector into its appropriate 2D form (rows are v0, columns v1)
        matmarg = sum(mat,dims=1) #Marginalize over one dimension to yield P(v0) = Σ_v1 P(v0,v1)
        U2[:,i] .= matmarg[1,:] #Save the 1D vector P(v0) as a column of our matrix
    end

    # Use a QR decomposition to orthogonalize the columns of the linear map
    U = Matrix(LinearAlgebra.qr(U2).Q)
    #Transpose U to yield a matrix whose rows are our basis vectors
    U = broadcast(real,transpose(U))

   return U

end

#############################################
# TRANSFORM YOUR MPO IN THE NEW BASIS: ######
#############################################

function TransformH(H,sitesNewBasis,U)

    # Arguments: 
        # H -> MPO master equation in our original occupation basis
        # sitesNewBasis -> the ITensor site indices for our MPO in our new basis
        # U -> a matrix whose rows are our new basis vectors
    # Returns:
	    # H2 -> MPO master equation in our new transformed basis

    # Create your ITensor MPO, with the appropriate site indices:
    H2 = MPO(sitesNewBasis)
    sitesOccBasis = firstsiteinds(H)

    ## Transform each site of the MPO:
    for i=1:(L+1) ## Loop over the voltage at each node, v_i
        Hten = H[i] #Start with the corresponding site of the MPO in the original basis
        Hten *= ITensor(U,sitesNewBasis[i],sitesOccBasis[i]) #Use the linear map to transform the incoming index
        Hten *= prime(ITensor(U,sitesNewBasis[i],sitesOccBasis[i]),plev=0,1) #Use the linear map to transform the outgoing index
        H2[i] = Hten #Populate the correct site of your new master equation
    end

    return H2

end

#############################################
# TRANSFORM YOUR MPS IN THE NEW BASIS: ######
#############################################

function TransformPsi(Psi,sitesNewBasis,U)

    # Arguments: 
        # Psi -> MPS in our original occupation basis
        # sitesNewBasis -> the ITensor site indices for our MPS in our new basis
        # U -> a matrix whose rows are our new basis vectors
    # Returns:
	    # Psi2 -> MPS state vector in our new transformed basis

    # Create your new ITensor MPS, with the appropriate site indices:
    Psi2 = MPS(sitesNewBasis)
    sitesOccBasis = siteinds(Psi)

    ## Transform each site of the MPO:
    for i=1:(L+1) ## Loop over the voltage at each node, v_i
        PsiTen = Psi[i] #Start with the corresponding site of the MPO in the original basis
        PsiTen *= ITensor(U,sitesNewBasis[i],sitesOccBasis[i]) #Use the linear map to transform the basis
        Psi2[i] = PsiTen #Populate the correct site of your new master equation
    end

    return Psi2

end

#########################################################
# TRANSFORM YOUR MPS BACK TO THE OCCUPATION BASIS: ######
#########################################################

function TransformPsiBack(Psi2,sitesOccBasis,U)

    # Arguments: 
        # Psi2 -> MPS in our transformed basis
        # sitesOccBasis -> the ITensor site indices for our MPS in our occupation basis
        # U -> a matrix whose rows are our new basis vectors
    # Returns:
	    # Psi -> MPS state vector in the original occupation basis

    # Create your new ITensor MPS, with the appropriate site indices:
    Psi = MPS(sitesOccBasis)
    sitesNewBasis = siteinds(Psi2)

    ## Transform each site of the MPO:
    for i=1:(L+1) ## Loop over the voltage at each node, v_i
        PsiTen = Psi2[i] #Start with the corresponding site of the MPO in the original basis
        PsiTen *= ITensor(U,sitesNewBasis[i],sitesOccBasis[i]) #Use the linear map to transform the basis
        Psi[i] = PsiTen #Populate the correct site of your new master equation
    end

    return Psi

end
