module InviscidBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Utils
using ..Poly

export Inviscid, lmn_t, lmn_p, s, t

struct Inviscid; end

Inviscid(N; kwargs...) = Basis{Inviscid}(;N, BC=InviscidBC(), kwargs...)

function t(::Basis{Inviscid}, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

function s(::Basis{Inviscid}, l,m,n,r) 
    fac = sqrt(5+2l+4n)/sqrt(4l*(l+1)*(n+1)^2)
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end

function lmn_p(b::Basis{Inviscid})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in 0:((N-l+1)÷2-1) if abs(m)<=l] 
    end
end

function lmn_t(b::Basis{Inviscid})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:N for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:N for m in ms for n in 0:((N-l)÷2) if abs(m)<=l] 
    end
end

lpmax(b::Basis{Inviscid}) = b.N-1
ltmax(b::Basis{Inviscid}) = b.N



function lmn_upol(N, ms = -N:N, ns = false) 
    if ns != false
        [(l,m,n) for l in 1:(N-1) for m in ms for n in ns if abs(m)<=l]
    else
        [(l,m,n) for l in 1:(N-1) for m in ms for n in 0:((N-l+1)÷2-1) if abs(m)<=l] 
    end
end

function lmn_utor(N, ms = -N:N, ns = false) 
    if ns != false
        [(l,m,n) for l in 1:N for m in ms for n in ns if abs(m)<=l]
    else
        [(l,m,n) for l in 1:N for m in ms for n in 0:((N-l)÷2) if abs(m)<=l] 
    end
end

n(N) = (2N^3+9N^2+7N)÷6
np(N) = ((-1)^N*(3 + (-1)^N*(-3+2N*(-1+N*(3+N)))))÷12
nt(N) = ((-1)^N*(-3 + (-1)^N*(3 + 2*N*(2 + N)*(4 + N))))÷12

function nlp(N,l)
    ((-1)^N*(3 - 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(1 - 6*l - 4*l^2 + 6*(2 + l)*N))))÷12
end

function nlt(N,l)
    ((-1)^N*(-3 + 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(13 - 4*l^2 + 6*(2 + l)*N))))÷12
end

lmn2k_p(l,m,n,N) = nlp(N,l-1) + (l+m)*((N-l+1)÷2) + n + 1
lmn2k_t(l,m,n,N) = nlt(N,l-1) + (l+m)*((N-l)÷2+1) + n + 1

function lmn_upol_l(N, ms = -N:N, ns=0)
    lmn = lmn_upol(N,ms,ns)
    lmnk = Vector{NTuple{4,Int}}[]
    L = N-1
    for _ in 1:L
        push!(lmnk,NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l,m,n = lmn[k]
        push!(lmnk[l], (k,l,m,n))
    end
    return lmnk
end

function lmn_utor_l(N, ms = -N:N, ns=0)
    lmn = lmn_utor(N,ms,ns)
    lmnk = Vector{NTuple{4,Int}}[]
    L = N
    for _ in 1:L
        push!(lmnk,NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l,m,n = lmn[k]
        push!(lmnk[l], (k,l,m,n))
    end
    return lmnk
end

_coriolis_tt(l,m; Ω = 2) = Ω*im*m/(l*(l+1))
_coriolis_ss(l,m; Ω = 2) = _coriolis_tt(l,m; Ω)

function _coriolis_ts(l,l2,m,m2,n,n2; Ω = 2.0) 
    if (m != m2)
        return nothing
    end
    if (l==l2+1) && (n==n2)
        
        return -Ω*sqrt((l^2-1)/(4l^2-1))*sqrt((l-m)*(l+m))/l
    elseif (l==l2-1) && (n==n2+1)
        return -Ω*sqrt(l*(l-m+1)*(l+m+1)/((2+l)*(2l+1)*(2l+3)))*(l+2)/(l+1)
    end
    return nothing
end

function _coriolis_st(l,l2,m,m2,n,n2; Ω = 2.0)
    aij = _coriolis_ts(l2,l,m2,m,n2,n; Ω)
    if !isnothing(aij)
        aij = -aij
    end
    return aij
end

@inline function _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0, thresh=sqrt(eps())) where T

    is,js,aijs = Int[], Int[], Complex{T}[]

    for ilmn in lmn_p_l
        for (i,l,m,n) in ilmn
            aij = _coriolis_ss(T(l),T(m); Ω)
            appendit!(is,js,aijs,i,i,aij; thresh)
            lrange = (l== 1) ? [2] : ((l==(N)) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_t_l,lrange)
                for (j,l2,m2,n2) in lmn2
                    aij = _coriolis_st(T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                    appendit!(is,js,aijs, i,j+np, aij)
                end
            end
        end
    end

    return is, js, aijs
end
@inline function _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0 ) where T

    is,js,aijs = Int[], Int[], Complex{T}[]

    for ilmn in lmn_t_l
        for (i,l,m,n) in ilmn
            aij = _coriolis_tt(T(l),T(m); Ω)
            appendit!(is,js,aijs, i+np,i+np, aij)

            lrange = (l== 1) ? [2] : ((l>=N-1) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_p_l,lrange)
                for (j,l2,m2,n2) in lmn2
                   aij =  _coriolis_ts(T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                   appendit!(is,js,aijs, i+np,j, aij)
                end
            end
        end
    end

    return is, js, aijs
end

function rhs_coriolis(N,m; ns = false, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m,ns)
    lmn_t = lmn_utor(N,m,ns)
    lmn_p_l = lmn_upol_l(N,m,ns)
    lmn_t_l = lmn_utor_l(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is, js, aijs = _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω )
    is2, js2, aijs2 = _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω )

    append!(is,is2)
    append!(js,js2)
    append!(aijs,aijs2)

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

function lhs(N,m; ns = false, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m, ns)
    lmn_t = lmn_utor(N,m, ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt


    LHS = sparse(I(nu)*one(Complex{T}))
    return LHS

end

end