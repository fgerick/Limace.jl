module InsulatingMFBasis

using SparseArrays

#lmns

function lmn_bpol(N, ms = 0:(N-1), ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:(N-1) for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:(N-1) for n in 1:(N-l+1)÷2 if abs(m)<=l] 
    end
end

function lmn_btor(N, ms = 0:(N-2), ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:(N-2) for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:(N-2) for n in 1:((N-l)÷2)  if abs(m)<=l]  
    end
end

#inner products

@inline function _inner_tt(is,js,aijs, i,j, l,n,n2)
    if n==n2 
        aij = one(l)
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    elseif (n==n2+1)
        # n = n-1
        aij = -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n)) + 3/(5 + 2*l + 4*(-1 + n)))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
        #symmetric
        push!(is,j)
        push!(js,i)
        push!(aijs,aij)

    end
        
    return nothing
end

@inline function _inner_ss(is,js,aijs, i,j, l,n,n2)
    
   
    if n==n2 
        aij = one(l)
        #also add the all space contribution here (only nonzero for n==n2==1):
        if n==1
            aij += l^2*(1 + l)*(5 + 2l)^2
        end
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    elseif (n==n2+1) 
        if n2 == 1
            aij = -(((1 + 2l)*(9 + 2l))/(sqrt(2*(7 + 2l)*(9 + 2l)*(6 + l*(11 + 6l)))))
            
        else
            aij = -((5 + 2l + 4n)*(9 + 2l + 4n))/2sqrt((-3 + 2l + 4n)*(-1 + 2l + 4n)*(3 + 2l + 4n)*(5 + 2l + 4n))
        end
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
        #symmetric
        push!(is,j)
        push!(js,i)
        push!(aijs,aij)
    end
        
    return nothing
end


function lhs(N,m; ns = 0, Ω::T = 1.0) where T
    lmn_p = lmn_bpol(N,m,ns)
    lmn_t = lmn_btor(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        _inner_ss(is,js,aijs,i,i, T(l),T(n),T(n))
        if n>1
            _inner_ss(is,js,aijs,i,i-1, T(l),T(n),T(n-1))
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        _inner_tt(is,js,aijs, i+np,i+np, T(l),T(n),T(n))
        if n>1
            _inner_tt(is,js,aijs, i+np,i-1+np, T(l),T(n),T(n-1))
        end
    end

    LHS = sparse(is,js,aijs, nu, nu)
    return LHS

end

#diffusion

@inline function _diffusion_tt(is,js,aijs, i,j, l,m,n,n2; η = 1.0) 
    
   
    if n==n2 
        aij = -((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,η*aij)
    end
        
    return nothing
end


@inline function _diffusion_ss(is,js,aijs, i,j, l,m,n,n2; η = 1.0) 
    
   
    if n==n2 
        if n==1
            aij = -(((3 + 2*l)*(5 + 12*l + 4*l^2)^2)/((5 + 2*l)*(6 + l*(11 + 6*l)))) #-((3 + 4*l*(2 + l))^2/(6 + l*(11 + 6*l)))
        else
            aij = -((-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/2
        end
        push!(is,i)
        push!(js,j)
        push!(aijs,η*aij)
    end

    return nothing
end

function rhs_diffusion(N,m; ns = 0, η::T = 1.0) where T
    lmn_p = lmn_bpol(N,m,ns)
    lmn_t = lmn_btor(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = T.((l,m,n))
        _diffusion_ss( is,js,aijs,i,i, l,m,n,n; η)
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = T.((l,m,n))
        _diffusion_tt( is,js,aijs, i+np,i+np, l,m,n,n; η)  
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end



end