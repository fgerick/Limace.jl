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
    params=Dict(:σw => σw, :σf => σf, :h => h)
    return Basis{ThinWall}(;N, V=Sphere(), BC=NoBC(), params,  kwargs...)
end

s(::Type{Basis{ThinWall}}, V::Volume, l,m,n,r) = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{ThinWall}}, V::Volume, l,m,n,r) = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{ThinWall},l) = 1:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{ThinWall},l) = 1:((b.N-l)÷2+1)

#10.1103/PhysRevE.88.053010
@inline function bcs_p(b::Basis{ThinWall}) 
    @inline _s = (l,n,r) -> r*s(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (
          @inline((l,n) -> σw*h/σf*(∂(r->∂(r->r*_s(l,n,r),r), r1) - l*(l+1)/r1^2*_s(l,n,r1)) + _s(l,n,r1)*l/r1 + ∂(r->r*_s(l,n,r),r1)*(1 + l*h/r1)), 
          )
    return fs
end

#10.1103/PhysRevE.88.053010
@inline function bcs_t(b::Basis{ThinWall}) 
    @inline _t = (l,n,r) -> r*t(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (@inline((l,n) -> σw/σf*h*∂(r->r*_t(l,n,r),r1) + _t(l,n,r1)), )
    return fs
end


lpmax(b::Basis{ThinWall}) = b.N
ltmax(b::Basis{ThinWall}) = b.N


end