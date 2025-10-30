module Bases

using DocStringExtensions

export BoundaryCondition, NoBC, InviscidBC, NoSlipBC, PerfectlyConductingBC, InsulatingBC, DirichletBC
export Volume, Sphere, SphericalShell
export LimaceBasis, Basis, BasisElement, isaxisymmetric, Helmholtz, Poloidal, Toroidal

import Base: length

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
abstract type BoundaryCondition end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct NoBC <: BoundaryCondition end
"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct InviscidBC <: BoundaryCondition end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct NoSlipBC <: BoundaryCondition end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct PerfectlyConductingBC <: BoundaryCondition end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct InsulatingBC <: BoundaryCondition end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
struct DirichletBC <: BoundaryCondition end

abstract type Volume end

Base.@kwdef struct Sphere <: Volume
    r1::Float64 = 1.0
end

struct SphericalShell <: Volume
    r0::Float64
    r1::Float64
end

import Base: getproperty

function Base.getproperty(s::Sphere, sym::Symbol)
    if sym === :r0
        return 0.0
    else
        return getfield(s, sym)
    end
end

abstract type LimaceBasis end

"""
$(TYPEDEF)

- `N::Int`: truncation degree
- `m::UnitRange{Int}`: spherical harmonic orders, default `-N:N``
- `n::UnitRange{Int}`: radial degrees, default `0:0` to make `n = n(N,l)`.
- `BC::BoundaryCondition`: boundary condition, default `NoBC()`
- `V::Vol`: volume, default `Sphere()`
- `params::Dict{Symbol,Float64}`: additional parameters, default empty dictionary

"""
Base.@kwdef struct Basis{T,Vol<:Volume} <: LimaceBasis
    N::Int #truncation degree
    m::UnitRange{Int} = -N:N #spherical harmonic orders
    n::UnitRange{Int} = 0:0 #radial degrees, default 0:0 to make n = n(N,l).
    BC::BoundaryCondition = NoBC()
    V::Vol = Sphere()
    params::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
end

Basis{T,V}(N::Int, m::Int, args...) where {T,V<:Volume} = Basis{T,V}(N, m:m, args...)


isaxisymmetric(b::Basis) = length(b.m) == 1


@inline function length(b::Basis)
    lmn_ps = lmn_p(b)
    lmn_ts = lmn_t(b)

    np = length(lmn_ps)
    nt = length(lmn_ts)
    nu = np + nt
    return nu
end

@inline function np(b::Basis)
    return length(lmn_p(b))
end

@inline function nt(b::Basis)
    return length(lmn_t(b))
end

@inline function nrange_p(b::Basis, l)
    if b.n == 0:0
        return _nrange_p(b, l)
    else
        return b.n
    end
end

@inline function nrange_t(b::Basis, l)
    if b.n == 0:0
        return _nrange_t(b, l)
    else
        return b.n
    end
end

@inline function _nrange_p(b::Basis, l)
    @error "define _nrange_p(b::Basis, l)!"
end

@inline function _nrange_t(b::Basis, l)
    @error "define _nrange_t(b::Basis, l)!"
end

@inline function nrange_p_bc(b::T, l) where {T<:Basis}
    nrange = nrange_p(b, l)
    if typeof(b.BC) != NoBC
        return nrange
    else
        nbc = length(bcs_p(b))
        return first(nrange):(last(nrange)-nbc)
    end
end

@inline function nrange_t_bc(b::T, l) where {T<:Basis}
    nrange = nrange_t(b, l)
    if typeof(b.BC) != NoBC
        return nrange
    else
        nbc = length(bcs_t(b))
        return first(nrange):(last(nrange)-nbc)
    end
end

@inline function lpmax(b::Basis)
    @error "define"
end

@inline function ltmax(b::Basis)
    @error "define"
end

function _lmn_l(lmn, L::Int)
    lmnk = Vector{NTuple{4,Int}}[]
    for _ in 1:L
        push!(lmnk, NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l, m, n = lmn[k]
        push!(lmnk[l], (k, l, m, n))
    end
    return lmnk
end

function lmn_t_l(b::Basis)
    lmn = lmn_t(b)
    L = ltmax(b)
    return _lmn_l(lmn, L)
end

function lmn_p_l(b::Basis)
    lmn = lmn_p(b)
    L = lpmax(b)
    return _lmn_l(lmn, L)
end

function lmn2k_dict(lmns)
    return Dict(lmn => i for (i, lmn) in enumerate(lmns))
end

lmn2k_p_dict(b::Basis) = lmn2k_dict(lmn_p(b))
lmn2k_t_dict(b::Basis) = lmn2k_dict(lmn_t(b))

function _lmn2cdeg_p(b::Basis, l, m, n)
    for N in 1:b.N
        _m = length(b.m) == 1 ? b.m : -N:N
        _n = b.n
        if (l,m,n) ∈ lmn_p(typeof(b)(;N, m=_m, n=_n))
            return N
        end
    end
    return nothing
end

function _lmn2cdeg_t(b::Basis, l, m, n)
    for N in 1:b.N
        _m = length(b.m) == 1 ? b.m : -N:N
        _n = b.n
        if (l,m,n) ∈ lmn_t(typeof(b)(;N, m=_m, n=_n))
            return N
        end
    end
    return nothing
end

function t(::Type{Basis}, V::Volume, l, m, n, r)
end

function s(::Type{Basis}, V::Volume, l, m, n, r)
end

t(b::T, l, m, n, r) where {T<:Basis} = t(T, b.V, l, m, n, r)
s(b::T, l, m, n, r) where {T<:Basis} = s(T, b.V, l, m, n, r)


function bcs_p(b::Basis)
    @error "implement"
end

function bcs_t(b::Basis)
    @error "implement"
end


function lmn_p(b::Basis)
    N, ms, ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l, m, n) for l in 1:lpmax(b) for m in ms for n in ns if abs(m) <= l]
    else
        return [(l, m, n) for l in 1:lpmax(b) for m in ms for n in nrange_p(b, l) if abs(m) <= l]
    end
end

function lmn_t(b::Basis)
    N, ms, ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l, m, n) for l in 1:ltmax(b) for m in ms for n in ns if abs(m) <= l]
    else
        return [(l, m, n) for l in 1:ltmax(b) for m in ms for n in nrange_t(b, l) if abs(m) <= l]
    end
end

function lmn_p_bc(b::Basis)
    N, ms, ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l, m, n) for l in 1:lpmax(b) for m in ms for n in ns if abs(m) <= l]
    else
        return [(l, m, n) for l in 1:lpmax(b) for m in ms for n in nrange_p_bc(b, l) if abs(m) <= l]
    end
end

function lmn_t_bc(b::Basis)
    N, ms, ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l, m, n) for l in 1:ltmax(b) for m in ms for n in ns if abs(m) <= l]
    else
        return [(l, m, n) for l in 1:ltmax(b) for m in ms for n in nrange_t_bc(b, l) if abs(m) <= l]
    end
end

abstract type Helmholtz end
struct Poloidal <: Helmholtz end
struct Toroidal <: Helmholtz end


"""
$(TYPEDEF)

- `TB<:Basis`: basis type
- `PT<:Helmholtz`: Helmholtz type, either `Poloidal` or `Toroidal`
- `lmn::NTuple{3,Int}`: spherical harmonic degree, order and radial degree, i.e `(l, m, n)`
- `factor::T`: factor, default `1.0`, can be used to scale the basis element

## Example
```julia
b = Insulating(10)
B0 = BasisElement(b, Poloidal, (1, 0, 0))
B1 = BasisElement(b, Toroidal, (2, 1, 2), 2.0)
```

Or without constructing a `Basis` object:
```julia
l, m, n = 1, 0, 1
B0 = BasisElement(Basis{Insulating,Sphere}, Poloidal, (l, m, n))
```

"""
struct BasisElement{TB<:Basis,PT<:Helmholtz,T<:Number}
    lmn::NTuple{3,Int}
    factor::T
end

import Base: length, iterate

length(u::BasisElement) = 1
iterate(u::BasisElement) = (u, nothing)
iterate(u::BasisElement, ::Any) = nothing

BasisElement(::TB, ::Type{PT}, lmn::NTuple{3,Int}, factor::T=1.0) where {TB<:Basis,PT<:Helmholtz,T<:Number} = BasisElement{TB,PT,T}(lmn, factor)
BasisElement(::Type{TB}, ::Type{PT}, lmn::NTuple{3,Int}, factor::T=1.0) where {TB<:Basis,PT<:Helmholtz,T<:Number} = BasisElement{TB,PT,T}(lmn, factor)


s(b::BasisElement{T,Poloidal}, V::Volume, r) where {T} = s(T, V, b.lmn..., r)
t(b::BasisElement{T,Toroidal}, V::Volume, r) where {T} = t(T, V, b.lmn..., r)


import Base: +, -, *

-(u::BasisElement{TB,PT}) where {TB<:Basis,PT<:Helmholtz} = BasisElement(TB, PT, u.lmn, -u.factor)
*(x::Number, u::BasisElement{TB,PT}) where {TB<:Basis,PT<:Helmholtz} = BasisElement(TB, PT, u.lmn, x * u.factor)

end #module
