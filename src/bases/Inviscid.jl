module InviscidBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, _lmn2cdeg_p, _lmn2cdeg_t
import ..Limace: inertial, _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!

export Inviscid
export PerfectlyConducting #same conditions on t and s

struct Inviscid; end
const PerfectlyConducting = Inviscid

Inviscid(N; kwargs...) = Basis{Inviscid}(;N, BC=InviscidBC(), V=Sphere(), kwargs...)

"""
$(TYPEDSIGNATURES)

```math
t_{l,n,m}(r) = f_{l,n}r^l J_n^{(0,l+1/2)}(2r^2-1)
```
with 

```math
f_{l,n} = \\sqrt{\\frac{3+2l+4n}{l(l+1)}}
```

[livermore_compendium_2014](@citet) (5.1), normalized to unit energy.
"""
@inline function t(::Type{Basis{Inviscid}}, V::Volume, l,m,n,r)
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

"""
$(TYPEDSIGNATURES)

```math
s_{l,n,m}(r) = f_{l,n}(1-r^2)r^l J_n^{(1,l+1/2)}(2r^2-1)
```

with
```math
f_{l,n} = \\sqrt{\\frac{5+2l+4n}{4l(l+1)(n+1)^2}}
```

[livermore_compendium_2014](@citet) (5.6), normalized to unit energy. 
"""
@inline function s(::Type{Basis{Inviscid}}, V::Volume, l,m,n,r)
    fac = sqrt(5+2l+4n)/sqrt(4l*(l+1)*(n+1)^2)
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end

@inline _nrange_p(b::Basis{Inviscid},l) = 0:((b.N-l+1)÷2-1)
@inline _nrange_t(b::Basis{Inviscid},l) = 0:((b.N-l)÷2)

@inline lpmax(b::Basis{Inviscid}) = b.N
@inline ltmax(b::Basis{Inviscid}) = b.N

_lmn2cdeg_p(u::Basis{Inviscid}, l,m,n) = l+2n+1
_lmn2cdeg_t(u::Basis{Inviscid}, l,m,n) = l+2n


n(N) = (2N^3+9N^2+7N)÷6

_np(N::Int) = ((-1)^N*(3 + (-1)^N*(-3+2N*(-1+N*(3+N)))))÷12

@inline function _np(N::Int, m::Int) 
    if (m == 0)
        m+=1
    end
    return (N-abs(m)+1)^2 ÷ 4
end

@inline function np(b::Basis{Inviscid})
    if isaxisymmetric(b)
        return _np(b.N,first(b.m))
    else
        return _np(b.N)
    end
end

_nt(N) = ((-1)^N*(-3 + (-1)^N*(3 + 2*N*(2 + N)*(4 + N))))÷12

@inline function _nt(N::Int, m::Int)
    if (m==0)
        m+=1
    end
    return (N-abs(m)+2)^2 ÷ 4
end

@inline function nt(b::Basis{Inviscid})
    if isaxisymmetric(b)
        return _nt(b.N,first(b.m))
    else
        return _nt(b.N)
    end
end


function nlp(N,l)
    ((-1)^N*(3 - 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(1 - 6*l - 4*l^2 + 6*(2 + l)*N))))÷12
end

function nlt(N,l)
    ((-1)^N*(-3 + 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(13 - 4*l^2 + 6*(2 + l)*N))))÷12
end

lmn2k_p(l,m,n,N) = nlp(N,l-1) + (l+m)*((N-l+1)÷2) + n + 1
lmn2k_t(l,m,n,N) = nlt(N,l-1) + (l+m)*((N-l)÷2+1) + n + 1


### Coriolis force analytically:


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

function _coriolis_poloidal_poloidal!(b::Basis{Inviscid}, is, js, aijs, lmn2k_p, l, m, r, wr, Ω::T) where T
    for n in nrange_p(b,l)
        aij = _coriolis_ss(T(l), T(m); Ω)
        appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n)], aij)
    end
    return nothing
end

function _coriolis_poloidal_toroidal!(b::Basis{Inviscid}, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω::T) where T
    for n in nrange_p(b,l), n2 in nrange_t(b,l2)
        aij = _coriolis_st(T(l),T(l2),T(m),T(m),T(n),T(n2); Ω)
        appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_t[(l2,m,n2)] + _np, aij)
    end
    return nothing
end


function _coriolis_toroidal_toroidal!(b::Basis{Inviscid}, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω::T) where T
    for n in nrange_t(b,l)
        aij = _coriolis_tt(T(l), T(m); Ω)
        appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n)] + _np, aij)
    end
    return nothing
end

function _coriolis_toroidal_poloidal!(b::Basis{Inviscid}, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω::T) where T
    for n in nrange_t(b,l), n2 in nrange_p(b,l2)
        aij = _coriolis_ts(T(l),T(l2),T(m),T(m),T(n),T(n2); Ω)
        appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_p[(l2,m,n2)], aij)
    end
    return nothing
end


function inertial(b::Basis{Inviscid})
    return one(typeof(b.V.r1))*I
end

end