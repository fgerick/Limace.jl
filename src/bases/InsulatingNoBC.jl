module InsulatingBasisNoBC

using SparseArrays


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis, ..InviscidBasis, ..InsulatingBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t
using ..Poly: ∂

export InsulatingNoBC

struct InsulatingNoBC end

InsulatingNoBC(N; kwargs...) = Basis{InsulatingNoBC}(; N, BC=NoBC(), V=Sphere(), kwargs...)

s(::Type{Basis{InsulatingNoBC}}, V::Volume, l,m,n,r)  = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{InsulatingNoBC}}, V::Volume, l,m,n,r)  = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{InsulatingNoBC}, l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{InsulatingNoBC}, l) = 0:((b.N-l)÷2)

@inline lpmax(b::Basis{InsulatingNoBC}) = b.N
@inline ltmax(b::Basis{InsulatingNoBC}) = b.N

@inline function bcs_t(b::Basis{InsulatingNoBC})
    fs = (@inline((l, n) -> t(Basis{InsulatingNoBC}, b.V, l, 0, n, 1.0)),)
    return fs
end

@inline function bcs_p(b::Basis{InsulatingNoBC})
    fs = (@inline((l, n) -> ∂(r -> s(Basis{InsulatingNoBC}, b.V, l, 0, n, r), 1.0) + (l + 1) * s(Basis{InsulatingNoBC}, b.V, l, 0, n, 1.0)),)
    # fs = ()
    return fs
end

end