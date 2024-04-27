module InviscidBasis2

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s
import ..Limace: _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!

export Inviscid2

struct Inviscid2; end

Inviscid2(N; kwargs...) = Basis{Inviscid2}(;N, BC=InviscidBC(), kwargs...)

function t(::Type{Basis{Inviscid2}}, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

function s(::Type{Basis{Inviscid2}}, l,m,n,r)
    fac = sqrt(5+2l+4n)/sqrt(4l*(l+1)*(n+1)^2)
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end

_nrange_p(b::Basis{Inviscid2},l) = 0:((b.N-l+1)÷2-1)
_nrange_t(b::Basis{Inviscid2},l) = 0:((b.N-l)÷2)


function lmn_p(b::Basis{Inviscid2})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:(N-1) for m in ms for n in nrange_p(b,l) if abs(m)<=l] 
    end
end

function lmn_t(b::Basis{Inviscid2})
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:N for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:N for m in ms for n in nrange_t(b,l) if abs(m)<=l] 
    end
end

lpmax(b::Basis{Inviscid2}) = b.N-1
ltmax(b::Basis{Inviscid2}) = b.N


n(N) = (2N^3+9N^2+7N)÷6

_np(N::Int) = ((-1)^N*(3 + (-1)^N*(-3+2N*(-1+N*(3+N)))))÷12

function _np(N::Int, m::Int) 
    if (m == 0)
        m+=1
    end
    return (N-abs(m)+1)^2 ÷ 4
end

function np(b::Basis{Inviscid2})
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

function nt(b::Basis{Inviscid2})
    if isaxisymmetric(b)
        return _nt(b.N,first(b.m))
    else
        return _nt(b.N)
    end
end

end