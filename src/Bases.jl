module Bases

using DocStringExtensions

export BoundaryCondition, NoBC, InviscidBC, NoSlipBC, PerfectlyConductingBC, InsulatingBC
export Volume
export Basis, BasisElement, isaxisymmetric, Helmholtz, Poloidal, Toroidal
# export nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax

import Base: length

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
abstract type BoundaryCondition; end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct NoBC <: BoundaryCondition; end
"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct InviscidBC <: BoundaryCondition; end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct NoSlipBC <: BoundaryCondition; end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct PerfectlyConductingBC <: BoundaryCondition; end

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct InsulatingBC <: BoundaryCondition; end


Base.@kwdef struct Volume
    r0::Float64 = 0.0
    r1::Float64 = 1.0
end

Sphere(; r=1.0) = Volume(;r1=r)
SphericalShell(r0,r1) = Volume(; r0,r1)


"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
Base.@kwdef struct Basis{T}
	N::Int #truncation degree
    m::UnitRange{Int} = -N:N #spherical harmonic orders
	n::UnitRange{Int} = 0:0 #radial degrees, default 0:0 to make n = n(N,l).
	BC::BoundaryCondition = NoBC()
    V::Volume = Sphere()
    params::Dict{Symbol,Float64} = Dict{Symbol,Float64}()
end

Basis{T}(N::Int, m::Int, n::UnitRange{Int}, BC::BoundaryCondition, args...) where {T} = Basis{T}(N, m:m, n, BC, args...)

isaxisymmetric(b::Basis) = length(b.m) == 1


@inline function length(b::Basis)
	lmn_ps = lmn_p(b)
    lmn_ts = lmn_t(b)

    np = length(lmn_ps)
    nt = length(lmn_ts)
    nu = np+nt
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
        return _nrange_p(b,l)
    else
        return b.n
    end
end

@inline function nrange_t(b::Basis, l)
    if b.n == 0:0
        return _nrange_t(b,l)
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

@inline function nrange_p_bc(b::T, l) where T<:Basis
    nrange = nrange_p(b,l)
    if typeof(b.BC) != NoBC
        return nrange
    else
        nbc = length(bcs_p(b))
        return first(nrange):(last(nrange)-nbc)
    end
end

@inline function nrange_t_bc(b::T, l) where T<:Basis
    nrange = nrange_t(b,l)
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
        push!(lmnk,NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l,m,n = lmn[k]
        push!(lmnk[l], (k,l,m,n))
    end
    return lmnk
end

function lmn_t_l(b::Basis)
    lmn = lmn_t(b)
    L = ltmax(b)
    return _lmn_l(lmn,L)
end

function lmn_p_l(b::Basis)
    lmn = lmn_p(b)
    L = lpmax(b)
    return _lmn_l(lmn,L)
end

function lmn2k_dict(lmns)
	return Dict(lmn=>i for (i,lmn) in enumerate(lmns))
end

lmn2k_p_dict(b::Basis) = lmn2k_dict(lmn_p(b))
lmn2k_t_dict(b::Basis) = lmn2k_dict(lmn_t(b))

# function t(::Type{Basis}, l, m, n, r)
# end

# function s(::Type{Basis}, l, m, n, r)
# end

function t(::Type{Basis}, V::Volume, l, m, n, r)
end

function s(::Type{Basis}, V::Volume, l, m, n, r)
end

t(b::T, l, m, n, r) where T<:Basis = t(T, b.V, l, m, n, r)
s(b::T, l, m, n, r) where T<:Basis = s(T, b.V, l, m, n, r)


function bcs_p(b::Basis)
    @error "implement"
end

function bcs_t(b::Basis)
    @error "implement"
end


function lmn_p(b::Basis)
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:lpmax(b) for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:lpmax(b) for m in ms for n in nrange_p(b,l) if abs(m)<=l] 
    end
end

function lmn_t(b::Basis)
    N,ms,ns = b.N, b.m, b.n
    if ns != 0:0
        return [(l,m,n) for l in 1:ltmax(b) for m in ms for n in ns if abs(m)<=l]
    else
        return [(l,m,n) for l in 1:ltmax(b) for m in ms for n in nrange_t(b,l) if abs(m)<=l] 
    end
end

abstract type Helmholtz end
struct Poloidal <: Helmholtz; end
struct Toroidal <: Helmholtz; end


"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
struct BasisElement{TB<:Basis,PT<:Helmholtz,T<:Number}
    lmn::NTuple{3,Int}
    factor::T
end

BasisElement(::TB, ::Type{PT}, lmn::NTuple{3,Int}, factor::T) where {TB<:Basis, PT<:Helmholtz, T<:Number} = BasisElement{TB,PT,T}(lmn, factor)
BasisElement(::Type{TB}, ::Type{PT}, lmn::NTuple{3,Int}, factor::T) where {TB<:Basis, PT<:Helmholtz, T<:Number} = BasisElement{TB,PT,T}(lmn, factor)


# s(b::T, l, m, n, r) where T<:Basis = s(T,l,m,n,r)
# t(b::T, l, m, n, r) where T<:Basis = t(T,l,m,n,r)
s(b::BasisElement{T,Poloidal}, V::Volume, r) where T = s(T,V,b.lmn..., r)
t(b::BasisElement{T,Toroidal}, V::Volume, r) where T = t(T,V,b.lmn..., r)


end #module