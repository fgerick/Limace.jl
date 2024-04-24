module Bases

export BoundaryCondition, NoBC, InviscidBC, NoSlipBC, PerfectlyConductingBC, InsulatingBC
export Basis, isaxisymmetric, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
export nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t

abstract type BoundaryCondition; end

struct NoBC <: BoundaryCondition; end
struct InviscidBC <: BoundaryCondition; end
struct NoSlipBC <: BoundaryCondition; end
struct PerfectlyConductingBC <: BoundaryCondition; end
struct InsulatingBC <: BoundaryCondition; end

import Base: length

Base.@kwdef struct Basis{T}
	N::Int #truncation degree
    m::UnitRange{Int} = -N:N #spherical harmonic orders
	n::UnitRange{Int} = 0:0 #radial degrees, default 0:0 to make n = n(N,l).
	BC::BoundaryCondition = NoBC()
end

Basis{T}(N::Int, m::Int, n::UnitRange{Int}, BC::BoundaryCondition) where T = Basis{T}(N, m:m, n, BC)

# Basis{T}(N::Int; m::Int = 0, n::UnitRange{Int} = 0:0, BC::BoundaryCondition = NoBC()) where T = Basis{T}(N; m = m:m, n, BC)

isaxisymmetric(b::Basis) = length(b.m) == 1


function length(b::Basis)
	lmn_ps = lmn_p(b)
    lmn_ts = lmn_t(b)

    np = length(lmn_ps)
    nt = length(lmn_ts)
    nu = np+nt
	return nu
end

function np(b::Basis)
	return length(lmn_p(b))
end

function nt(b::Basis)
	return length(lmn_t(b))
end

@inline function nrange_p(b::Basis, l)
    if b.n == 0:0
        return _nrange_p(b,l)
    else
        return b.n
    end
	# @error "define nrange_p(b::Basis, l)!"
end

@inline function nrange_t(b::Basis, l)
    if b.n == 0:0
        return _nrange_t(b,l)
    else
        return b.n
    end
	# @error "define nrange_p(b::Basis, l)!"
end

function _nrange_p(b::Basis, l)
	@error "define _nrange_p(b::Basis, l)!"
end

function _nrange_t(b::Basis, l)
	@error "define _nrange_t(b::Basis, l)!"
end

function nrange_p_bc(b::Basis, l)
	@error "define nrange_p(b::Basis, l)!"
end

function nrange_t_bc(b::Basis, l)
	@error "define nrange_t(b::Basis, l)!"
end

# function lmn_p(b::T) where T<:Basis
# 	@error "lmn_p(b) not defined for basis of type $T"
# end

# function lmn_t(b::Basis)
# end


function lpmax(b::Basis)
end

function ltmax(b::Basis)
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

function t(::Type{Basis}, l, m, n, r)
end

function s(::Type{Basis}, l, m, n, r)
end

function bcs_p(::Type{Basis})
end

function bcs_t(::Type{Basis})
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



end #module