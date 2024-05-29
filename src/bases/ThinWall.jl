module ThinWallBC

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

export ThinWall

struct ThinWall; end

function ThinWall(N; σw=1.0, σf = 1.0, h = 1.0, kwargs...)
    params=(; σw, σf, h)
    return Basis{ThinWall, typeof(params)}(;N, V=Sphere(), BC=NoBC(), params,  kwargs...)
end

s(::Type{Basis{ThinWall,P}}, V::Volume, l,m,n,r) where P  = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{ThinWall,P}}, V::Volume, l,m,n,r) where P = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{ThinWall,P},l) where P = 0:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{ThinWall,P},l) where P = 0:((b.N-l)÷2+1)

@inline function bcs_p(b::Basis{ThinWall,P}) where P
    @inline _s = (l,n,r) -> r*s(Basis{ThinWall,P}, b.V, l, 0, n, r)
    (; r1) = b.V 
    (; h, σf, σw) = b.params
    fs = (
          @inline((l,n) -> σw*h/σf*(∂(r->∂(r->_s(l,n,r),r), r1) - l*(l+1)/r1^2*_s(l,n,r1)) + _s(l,n,r1)*l/r1 + ∂(r->_s(l,n,r),r1)*(1 + l*h/r1)), 
          )
    return fs
end

@inline function bcs_t(b::Basis{ThinWall,P}) where P
    @inline _t = (l,n,r) -> r*t(Basis{ThinWall,P}, b.V, l, 0, n, r)
    (; r1) = b.V 
    (; h, σf, σw) = b.params
    fs = (@inline((l,n) -> σw/σf*h*∂(r->_t(l,n,r),r1) + _t(l,n,r1)), )
    return fs
end


lpmax(b::Basis{ThinWall}) = b.N
ltmax(b::Basis{ThinWall}) = b.N


end