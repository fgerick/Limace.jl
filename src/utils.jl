module Utils

export BoundaryCondition, NoBC, InviscidBC, NoSlipBC, PerfectlyConductingBC, InsulatingBC
export Basis, isaxisymmetric, lmn_p_l, lmn_t_l, appendit!

abstract type BoundaryCondition; end

struct NoBC <: BoundaryCondition; end
struct InviscidBC <: BoundaryCondition; end
struct NoSlipBC <: BoundaryCondition; end
struct PerfectlyConductingBC <: BoundaryCondition; end
struct InsulatingBC <: BoundaryCondition; end

Base.@kwdef struct Basis{T}
	N::Int #truncation degree
    m::UnitRange{Int} = -N:N #spherical harmonic orders
	n::UnitRange{Int} = 0:0 #radial degrees, default 0:0 to make n = n(N,l).
	BC::BoundaryCondition = NoBC()
end

Basis{T}(N::Int, m::Int, n::UnitRange{Int}, BC::BoundaryCondition) where T = Basis{T}(N, m:m, n, BC)

# Basis{T}(N::Int; m::Int = 0, n::UnitRange{Int} = 0:0, BC::BoundaryCondition = NoBC()) where T = Basis{T}(N; m = m:m, n, BC)

isaxisymmetric(b::Basis) = length(b.m) == 1

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
    lmn = lmn_t(b)
    L = lpmax(b)
    return _lmn_l(lmn,L)
end


function appendit!(is, js, aijs, i, j, aij; thresh=sqrt(eps()))
    if !isnothing(aij) && (abs(aij) > thresh)
        push!(is, i)
        push!(js, j)
        push!(aijs, aij)
    end
end



end #module