module ViscousBasisNoBC

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis
using ..InviscidBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export ViscousNoBC

struct ViscousNoBC; end

ViscousNoBC(N; kwargs...) = Basis{ViscousNoBC, NamedTuple}(;N, V=Sphere(), BC=NoBC(), kwargs...)

s(::Type{Basis{ViscousNoBC, NamedTuple}}, V::Volume, l,m,n,r)  = s(Basis{Inviscid,NamedTuple}, V, l,m,n,r) 
t(::Type{Basis{ViscousNoBC, NamedTuple}}, V::Volume, l,m,n,r)  = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{ViscousNoBC, NamedTuple},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{ViscousNoBC, NamedTuple},l) = 0:((b.N-l)÷2)

@inline function bcs_p(b::Basis{ViscousNoBC, NamedTuple})
    fs = (
        #   @inline((l,n)->s(Basis{ViscousNoBC, NamedTuple}, b.V, l, 0, n, b.V.r1)),
          @inline((l,n)->∂(r->s(Basis{ViscousNoBC, NamedTuple}, b.V, l, 0, n, r), b.V.r1)), 
          )
    return fs
end

@inline function bcs_t(b::Basis{ViscousNoBC, NamedTuple})
    fs = (@inline((l,n)->t(Basis{ViscousNoBC, NamedTuple}, b.V, l, 0, n, b.V.r1)), )
    return fs
end


lpmax(b::Basis{ViscousNoBC, NamedTuple}) = b.N-1
ltmax(b::Basis{ViscousNoBC, NamedTuple}) = b.N


end