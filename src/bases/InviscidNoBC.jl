module InviscidBasisNoBC

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export InviscidNoBC

struct InviscidNoBC; end

InviscidNoBC(N; kwargs...) = Basis{InviscidNoBC}(;N, V=Sphere(), BC=NoBC(), kwargs...)

s(::Type{Basis{InviscidNoBC}}, V::Volume, l,m,n,r)  = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{InviscidNoBC}}, V::Volume, l,m,n,r)  = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{InviscidNoBC},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{InviscidNoBC},l) = 0:((b.N-l)÷2)

@inline function bcs_p(b::Basis{InviscidNoBC})
    fs = (@inline((l,n)->s(Basis{InviscidNoBC}, b.V, l, 0, n, 1.0)), )
    # return ()
    return fs
end

@inline function bcs_t(b::Basis{InviscidNoBC})
    return ()
end


lpmax(b::Basis{InviscidNoBC}) = b.N
ltmax(b::Basis{InviscidNoBC}) = b.N


end