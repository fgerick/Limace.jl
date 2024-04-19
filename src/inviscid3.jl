module InviscidBasis3

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Utils
using ..Poly

import ..Utils: lpmax, ltmax, lmn_t, lmn_p, nrange_p, nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t
import ..Limace: _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!

export Inviscid3

struct Inviscid3; end

Inviscid3(N; kwargs...) = Basis{Inviscid3}(;N, BC=NoBC(), kwargs...)

@inline function t(::Type{Basis{Inviscid3}}, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

@inline function s(::Type{Basis{Inviscid3}}, l,m,n,r)
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
    # return r^l*(jacobi(n+1,0,l+1/2, 2r^2-1) - jacobi(n,0,l+1/2, 2r^2-1))
end

nrange_p(b::Basis{Inviscid3},l) = 0:((b.N-l+1)÷2-1)
nrange_t(b::Basis{Inviscid3},l) = 0:((b.N-l)÷2)

nrange_p_bc(b::Basis{Inviscid3},l) = 0:((b.N-l+1)÷2-2)
nrange_t_bc(b::Basis{Inviscid3},l) = nrange_t(b, l)

@inline function bcs_p(::Type{Basis{Inviscid3}})
    fs = (@inline((l,n)->s(Basis{Inviscid3}, l, 0, n, 1.0)), )
    return fs
end

@inline function bcs_t(::Type{Basis{Inviscid3}})
    return ()
end

function lmn_p(b::Basis{Inviscid3})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in nrange_p(b,l) if abs(m)<=l] 
    end
end

function lmn_t(b::Basis{Inviscid3})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:N for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:N for m in ms for n in nrange_t(b,l) if abs(m)<=l] 
    end
end

lpmax(b::Basis{Inviscid3}) = b.N-1
ltmax(b::Basis{Inviscid3}) = b.N


n(N) = (2N^3+9N^2+7N)÷6

_np(N::Int) = ((-1)^N*(3 + (-1)^N*(-3+2N*(-1+N*(3+N)))))÷12

function _np(N::Int, m::Int) 
    if (m == 0)
        m+=1
    end
    return (N-abs(m)+1)^2 ÷ 4
end

function np(b::Basis{Inviscid3})
    if isaxisymmetric(b)
        return _np(b.N,first(b.m))
    else
        return _np(b.N)
    end
end

_nt(N) = ((-1)^N*(-3 + (-1)^N*(3 + 2*N*(2 + N)*(4 + N))))÷12

function _nt(N::Int, m::Int)
    if (m==0)
        m+=1
    end
    return (N-abs(m)+2)^2 ÷ 4
end

function nt(b::Basis{Inviscid3})
    if isaxisymmetric(b)
        return _nt(b.N,first(b.m))
    else
        return _nt(b.N)
    end
end

end