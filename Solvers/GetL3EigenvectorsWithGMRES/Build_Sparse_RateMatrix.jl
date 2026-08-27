using SparseArrays

#################################################################################
# Sparse rate-matrix builder, general in the chain length L.                    #
#                                                                                #
# This mirrors Solvers/GetEigenvectorsWithDMRG/BuildMPO.jl term-by-term: the "Series"   #
# branch here loops over the same L nodes, in the same order, applying the same #
# "b operators" / "t operators" gain/loss terms, using the SAME operator names  #
# (x+, x-, U+, U-, D+, D-, upperProj, lowerProj) as that file's AutoMPO calls.  #
# Read the two files side by side -- every `ampo .+=`/`ampo .-=` line in        #
# BuildMPO.jl has a matching `H += `/`H -= ` line here, on the same site        #
# indices. That correspondence is the point: this file is the independent      #
# sparse/exact solve that Figure 3 checks the tensor-network result against,    #
# so it should be checkable by eye against the MPO it's meant to reproduce.     #
#                                                                                #     #
#################################################################################

function site_op(op, site, nsites, ID)
    # Kron `op` into position `site` (1-indexed) among `nsites` identities.
    factors = Any[ID for _ in 1:nsites]
    factors[site] = op
    out = factors[1]
    for k=2:nsites
        out = kron(out, factors[k])
    end
    return out
end

function two_site_op(opA, siteA, opB, siteB, nsites, ID)
    # Kron opA into siteA and opB into siteB among nsites identities (siteA != siteB).
    factors = Any[ID for _ in 1:nsites]
    factors[siteA] = factors[siteA] * opA
    factors[siteB] = factors[siteB] * opB
    out = factors[1]
    for k=2:nsites
        out = kron(out, factors[k])
    end
    return out
end

function H_Mechanisms(L,n,Xplus,Xminus,Uplus,Uminus,Dplus,Dminus,upperProj,lowerProj,ID,kw)
    # Arguments:
        # L -> number of units in the chain (nsites = L+1)
        # n -> physical (per-site) dimension, 2M+1
        # Xplus,Xminus,Uplus,Uminus,Dplus,Dminus,upperProj,lowerProj,ID -> the
        #     per-site sparse operators built in L3_CMOS_System.jl, named to match
        #     Solvers/GetEigenvectorsWithDMRG/MPOoperators.jl
        # kw -> "Series" for a chain of L units connected in series (kw="SingleUnit"
        #     is just the L=1 case of "Series" and is not special-cased separately,
        #     matching how BuildMPO.jl's two branches are algebraically the same
        #     construction at L=1)
    # Returns:
        # H -> the sparse rate matrix (generator) for the chain, size n^(L+1) x n^(L+1)

    nsites = L + 1
    H = spzeros(n^nsites, n^nsites)

    if kw == "Series" || kw == "SingleUnit"

        for i=1:L #loop over the nodes i (1-indexed; node i and its gate node i+1)

            ########################################
            #### b operators (bottom NOT gate) #####
            ########################################

            ### GAIN TERMS ### (off diagonal)

            # λu+(vi-ve;v(i+1)):
            H += two_site_op(Xplus, i, Uplus, i+1, nsites, ID)
            # λd-(vi-ve;v(i+1)):
            H += two_site_op(Xplus*Dminus, i, Uminus, i+1, nsites, ID)
            # λu-(vi+ve;v(i+1)):
            H += two_site_op(Xminus, i, Uminus, i+1, nsites, ID)
            # λd+(vi+ve;v(i+1)):
            H += two_site_op(Xminus*Dplus, i, Uplus, i+1, nsites, ID)

            ### LOSS TERMS ### (diagonal)

            # λu+(vi;v(i+1)):
            H -= two_site_op(upperProj, i, Uplus, i+1, nsites, ID)
            # λd-(vi;v(i+1)):
            H -= two_site_op(Dminus*upperProj, i, Uminus, i+1, nsites, ID)
            # λu-(vi;v(i+1)):
            H -= two_site_op(lowerProj, i, Uminus, i+1, nsites, ID)
            # λd+(vi;v(i+1)):
            H -= two_site_op(Dplus*lowerProj, i, Uplus, i+1, nsites, ID)

            ########################################
            #### t operators (top NOT gate) ########
            ########################################

            ### GAIN TERMS ### (off diagonal)

            # λu+(v(i+1)-ve;vi):
            H += two_site_op(Uplus, i, Xplus, i+1, nsites, ID)
            # λd-(v(i+1)-ve;vi):
            H += two_site_op(Uminus, i, Xplus*Dminus, i+1, nsites, ID)
            # λu-(v(i+1)+ve;vi):
            H += two_site_op(Uminus, i, Xminus, i+1, nsites, ID)
            # λd+(v(i+1)+ve;vi):
            H += two_site_op(Uplus, i, Xminus*Dplus, i+1, nsites, ID)

            ### LOSS TERMS ### (diagonal)

            # λu+(v(i+1);vi):
            H -= two_site_op(Uplus, i, upperProj, i+1, nsites, ID)
            # λd-(v(i+1);vi):
            H -= two_site_op(Uminus, i, Dminus*upperProj, i+1, nsites, ID)
            # λu-(v(i+1);vi):
            H -= two_site_op(Uminus, i, lowerProj, i+1, nsites, ID)
            # λd+(v(i+1);vi):
            H -= two_site_op(Uplus, i, Dplus*lowerProj, i+1, nsites, ID)

        end

    end

    return H
end
