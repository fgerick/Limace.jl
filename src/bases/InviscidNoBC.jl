module InviscidBasisNoBC

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export InviscidNoBC

struct InviscidNoBC; end

InviscidNoBC(N; kwargs...) = Basis{InviscidNoBC}(;N, BC=NoBC(), kwargs...)

t(::Type{Basis{InviscidNoBC}}, l,m,n,r)  = t(Basis{Unconstrained}, l,m,n,r) 
s(::Type{Basis{InviscidNoBC}}, l,m,n,r)  = s(Basis{Unconstrained}, l,m,n,r) 

@inline _nrange_p(b::Basis{InviscidNoBC},l) = 1:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{InviscidNoBC},l) = 0:((b.N-l)÷2)

# @inline nrange_p_bc(b::Basis{InviscidNoBC},l) = 0:((b.N-l+1)÷2-1)
# @inline nrange_t_bc(b::Basis{InviscidNoBC},l) = nrange_t(b, l)

@inline function bcs_p(::Type{Basis{InviscidNoBC}})
    fs = (@inline((l,n)->s(Basis{InviscidNoBC}, l, 0, n, 1.0)), )
    # return ()
    return fs
end

@inline function bcs_t(::Type{Basis{InviscidNoBC}})
    return ()
end


lpmax(b::Basis{InviscidNoBC}) = b.N-1
ltmax(b::Basis{InviscidNoBC}) = b.N


end