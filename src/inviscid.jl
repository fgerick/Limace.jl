module InviscidBasis

using SparseArrays
using LinearAlgebra

# lmn_upol(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:(N-1) for n in 0:(N-l+1)÷2 if abs(m)<=l]
# lmn_utor(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:N for n in 0:((N-l)÷2) if abs(m)<=l]

function lmn_upol(N, ms = 0:N, ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:(N-1) for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:(N-1) for n in 0:((N-l+1)÷2) if abs(m)<=l] 
    end
end

function lmn_utor(N, ms = 0:N, ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:N for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:N for n in 0:((N-l)÷2) if abs(m)<=l] 
    end
end

function lmn_upol_l(N, ms = 0:N, ns=0)
    # lmn = [[(l,m,n) for m in ms for n in 1:((N-l+1)÷2+1) if abs(m)<=l] for l in 1:(N-1)] 
    if (ns != 0)
        lmn = [[(l,m,n) for m in ms for n in ns if abs(m)<=l] for l in 1:(N-1)]
    else
        lmn = [[(l,m,n) for m in ms for n in 0:((N-l+1)÷2) if abs(m)<=l] for l in 1:(N-1)]
    end
    lmnk = Vector{NTuple{4,Int}}[]
    k=1
    for l in eachindex(lmn)
        push!(lmnk,NTuple{4,Int}[])
        for lmn in lmn[l]
            push!(lmnk[l], (k,lmn...))
            k+=1
        end
    end
    return lmnk
end

function lmn_utor_l(N, ms = 0:N, ns=0)
    if (ns != 0)
        lmn = [[(l,m,n) for m in ms for n in ns if abs(m)<=l] for l in 1:N]
    else
        lmn = [[(l,m,n) for m in ms for n in 0:((N-l)÷2) if abs(m)<=l] for l in 1:N]
    end
    lmnk = Vector{NTuple{4,Int}}[]
    k=1
    for l in eachindex(lmn)
        push!(lmnk,NTuple{4,Int}[])
        for lmn in lmn[l]
            push!(lmnk[l], (k,lmn...))
            k+=1
        end
    end
    return lmnk
end


_coriolis_tt(l,m; Ω = 2) = Ω*im*m/(l*(l+1))
_coriolis_ss(l,m; Ω = 2) = _coriolis_tt(l,m; Ω)

function _coriolis_tt(is,js,aijs, i,j, l,m; Ω = 2.0)
    push!(is,i)
    push!(js,j)
    push!(aijs, _coriolis_tt(l,m; Ω))
end

function _coriolis_ss(is,js,aijs, i,j, l,m; Ω = 2.0)
    push!(is,i)
    push!(js,j)
    push!(aijs, _coriolis_ss(l,m; Ω))
end


function _coriolis_ts(is,js,aijs, i,j, l,l2,m,m2,n,n2; Ω = 2.0) 
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

function _coriolis_st(is,js,aijs, i,j, l2,l,m2,m,n2,n; Ω = 2.0) 
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


@inline function _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0 ) where T

    is,js,aijs = Int[], Int[], Complex{T}[]

    for (l,ilmn) in enumerate(lmn_p_l)
        for (i,_,m,n) in ilmn
            _coriolis_ss(is,js,aijs, i,i, T(l),T(m); Ω)
            # for (j,l2,m2,n2) in lmn_p_l[l]
            #     if (m==m2) && (n==n2)
            #     end
            # end
            lrange = (l== 1) ? [2] : ((l==(N-1)) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_t_l,lrange)
                # l2,m2,n2 = T.(lmn_p[j])
                for (j,l2,m2,n2) in lmn2
                    _coriolis_st(is,js,aijs, i,j+np, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                end
            end
        end
    end

    return is, js, aijs
end
@inline function _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0 ) where T

    # np = length(lmn_p)
    # nt = length(lmn_t)
    # nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (l,ilmn) in enumerate(lmn_t_l)
        for (i,_,m,n) in ilmn
        #     # l,m,n = T.((l,m,n))
        #     for (j,l2,m2,n2) in lmn_t_l[l]
        #         if (m==m2) && (n==n2)
        #         # l2,m2,n2 = T.(lmn_t[j])
        #         _coriolis_tt(is,js,aijs, i+np,j+np, T(l),T(m); Ω)
        #         end
            _coriolis_tt(is,js,aijs, i+np,i+np, T(l),T(m); Ω)
            # end

            # lrange = l == 1 ? [1,2] : (l-1:l+1)
            lrange = (l== 1) ? [2] : ((l>=N-2) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_p_l,lrange)
                # l2,m2,n2 = T.(lmn_p[j])
                for (j,l2,m2,n2) in lmn2
                    _coriolis_ts(is,js,aijs, i+np,j, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                end
            end
        end
    end

    return is, js, aijs
end
function rhs_coriolis_2(N,m; ns = 0, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m,ns)
    lmn_t = lmn_utor(N,m,ns)
    lmn_p_l = lmn_upol_l(N,m,ns)
    lmn_t_l = lmn_utor_l(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    # is, js, aijs = rhs_coriolis_1(eachindex(lmn_p), lmn_p, lmn_t; Ω )
    # is2, js2, aijs2 = rhs_coriolis_2(eachindex(lmn_t), lmn_p, lmn_t; Ω )
    is, js, aijs = _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω )
    is2, js2, aijs2 = _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω )

    append!(is,is2)
    append!(js,js2)
    append!(aijs,aijs2)

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

function rhs_coriolis(N,m; ns = 0, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        # l,m,n = Float64.((l,m,n))
        push!(is,i)
        push!(js,i)
        push!(aijs,_coriolis_ss(T(l),T(m); Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            # l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_st(is,js,aijs, i,j+np, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        # l,m,n = Float64.((l,m,n))
        push!(is,i+np)
        push!(js,i+np)
        push!(aijs,_coriolis_tt(T(l),T(m); Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            # l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_ts(is,js,aijs, i+np,j, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
        end
    end

    RHS = sparse(is,js,aijs,nu,nu)
    return RHS

end

function lhs(N,m; ns = 0, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt


    LHS = sparse(I(nu)*one(Complex{T}))
    return LHS

end

end