lmn_upol_inviscid(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:(N-1) for n in 0:(N-l+1)÷2 if abs(m)<=l]
lmn_utor_inviscid(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:N for n in 0:((N-l)÷2) if abs(m)<=l]

# _coriolis_tt(l,m) = im*m/(l*(l+1))
# _coriolis_ss(l,m) = _coriolis_tt(l,m)

function _coriolis_ts_inviscid(is,js,aijs, i,j, l,l2,m,m2,n,n2; Ω = 2.0) 
    if (m != m2)
        return nothing
    end
    if (l==l2+1) && (n==n2)
        
        aij = -sqrt((l^2-1)/(4l^2-1))*sqrt((l-m)*(l+m))/l
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    elseif (l==l2-1) && (n==n2+1)
        aij = -sqrt(l*(l-m+1)*(l+m+1)/((2+l)*(2l+1)*(2l+3)))*(l+2)/(l+1)
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    end
    return nothing
end

function _coriolis_st_inviscid(is,js,aijs, i,j, l2,l,m2,m,n2,n; Ω = 2.0) 
    if (m != m2)
        return nothing
    end
    if (l==l2+1) && (n==n2)
        
        aij = sqrt((l^2-1)/(4l^2-1))*sqrt((l-m)*(l+m))/l
        push!(is,i)
        push!(js,j)
        push!(aijs, Ω*aij)
    elseif (l==l2-1) && (n==n2+1)
        aij = sqrt(l*(l-m+1)*(l+m+1)/((2+l)*(2l+1)*(2l+3)))*(l+2)/(l+1)
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    end
    return nothing
end

function rhs_inviscid(N,m; Ω = 2.0)
    lmn_p = lmn_upol_inviscid(N,m)
    lmn_t = lmn_utor_inviscid(N,m)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], ComplexF64[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        push!(is,i)
        push!(js,i)
        push!(aijs,_coriolis_ss(l,m; Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            _coriolis_st_inviscid(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        push!(is,i+np)
        push!(js,i+np)
        push!(aijs,_coriolis_tt(l,m; Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            _coriolis_ts_inviscid(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end

    RHS = sparse(is,js,aijs)
    return RHS

end