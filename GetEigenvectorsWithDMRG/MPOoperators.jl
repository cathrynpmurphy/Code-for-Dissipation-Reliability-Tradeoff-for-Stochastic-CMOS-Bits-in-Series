
## RAISING AND LOWERING OPERATORS ##

function ITensors.op!(Op::ITensor,
    ::OpName"x-", # SO CALLED "ac" because lowering operator
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M) # "for all a from 1 to the cap"
        Op[s'=>a,s=>a+1] = 1 # This constructs the operator which is a tensor of indices s x s', and where the upper diagonal is equal to the value (lowering operator)
    end
end

function ITensors.op!(Op::ITensor,
    ::OpName"x+", # SO CALLED "ac_dag" because raising operator
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M)
        Op[s'=>a+1,s=>a] =1 # This constructs the operator which is a tensor of indices s x s', and where the lower diagonal is equal to one (raising operator)
    end
end

## IDENTITY OPERATOR ##

function ITensors.op!(Op::ITensor,
    ::OpName"ID",
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M+1)
        Op[s'=>a,s=>a] = 1 # A diagonal of ones (identity matrix)
    end
end

## DIAGONAL RATE OPERATORS ##

function ITensors.op!(Op::ITensor,
    ::OpName"U+", # Uphill associated with PMOS, used to construct b^u_{+}, t^u_{+}, b^d_{+}, t^d_{+}
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M+1)
	mval = Mvals[a]
        Op[s'=>a,s=>a] = exp(vDD_vT/nval)*exp(-1*mval*ve_vT/nval) # This constructs the operator which is a tensor of indices s x s', and where the lower diagonal is equal to one (raising operator)
    end
end

function ITensors.op!(Op::ITensor,
    ::OpName"D+", # Downhill associated with PMOS, used to construct b^d_{+}, t^d_{+}
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M+1)
        mval = Mvals[a]
        Op[s'=>a,s=>a] = exp((-1/2)*ve_vT-vDD_vT)*exp(mval*ve_vT) # This constructs the operator which is a tensor of indices s x s', and where the lower diagonal is equal to one (raising operator)
    end
end

function ITensors.op!(Op::ITensor,
    ::OpName"U-", # Uphill associated with NMOS, used to construct b^u_{-}, t^u_{-}, b^d_{-}, t^d_{-}
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M+1)
	mval = Mvals[a]
        Op[s'=>a,s=>a] = exp(vDD_vT/nval)*exp(mval*ve_vT/nval) # This constructs the operator which is a tensor of indices s x s', and where the lower diagonal is equal to one (raising operator)
    end
end

function ITensors.op!(Op::ITensor,
    ::OpName"D-", # Downhill associated with NMOS, used to construct b^d_{+}, t^d_{+} 
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M+1)
        mval = Mvals[a]
        Op[s'=>a,s=>a] = exp((-1/2)*ve_vT-vDD_vT)*exp(-1*mval*ve_vT) # This constructs the operator which is a tensor of indices s x s', and where the lower diagonal is equal to one (raising operator)
    end
end

## OPERATORS FOR FINITE BOUNDARY CONDITIONS ##

function ITensors.op!(Op::ITensor,
    ::OpName"upperProj", 
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M)
        Op[s'=>a,s=>a] =1 # Zeros out the discrete voltage state at the positive boundary to keep probability from flowing out of bounds
    end
end

function ITensors.op!(Op::ITensor,
    ::OpName"lowerProj", 
    ::SiteType"Fock",
    s::Index)
    for a=1:(2*M)
        Op[s'=>a+1,s=>a+1] =1 # Zeros out the discrete voltage state at the negative boundary to keep probability from flowing out of bounds
    end
end
