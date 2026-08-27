function NormalizePss(Pss)

    # Arguments: 
        # Pss -> MPS of the steady state
    # Returns:
	    # PssNorm -> Steady State MPS normalized with the L1 norm <1|Pss> = 1

    sites = siteinds(Pss) #The ITensor site indices for the steady state MPS
    n = dim(sites[1]) #The physical dimension

    N1 = Pss[1]*ITensor(ones(1,n),sites[1]) #Marginalize: <1|Pss(v0)>
    for i=2:(L+1) #Loop over the sites of the MPS (the nodes i of the circuit)
        N1 *= Pss[i]*ITensor(ones(1,n),sites[i]) #Marginalize: <1|Pss(vi)>
    end

    Nval = scalar(N1)
    PssNorm = Pss ./ Nval

    return PssNorm,Nval
    
end
