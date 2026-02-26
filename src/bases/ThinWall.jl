module ThinWallBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis
using ..InviscidBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

import ..Limace: inertial, _inertial_ss, _inertial_tt
import ..Quadrature: rquad

export ThinWall

struct ThinWall; end

function ThinWall(N; σw=1.0, σf = 1.0, h = 0.0, μr = 1.0, kwargs...)
    params=Dict(:σw => σw, :σf => σf, :h => h, :μr => μr)
    return Basis{ThinWall,Sphere}(;N, V=Sphere(), BC=NoBC(), params,  kwargs...)
end

s(::Type{Basis{ThinWall,Sphere}}, V::Volume, l,m,n,r) = s(Basis{Unconstrained, Sphere}, V, l,m,n,r) 
t(::Type{Basis{ThinWall,Sphere}}, V::Volume, l,m,n,r) = t(Basis{Unconstrained, Sphere}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{ThinWall,Sphere},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{ThinWall,Sphere},l) = 0:((b.N-l)÷2)

#10.1103/PhysRevE.88.053010
@inline function bcs_p(b::Basis{ThinWall,Sphere}) 
    @inline _s = (l,n,r) -> r*s(Basis{ThinWall,Sphere}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw, μr = b.params[:h], b.params[:σf], b.params[:σw], b.params[:μr]
    fs = (
          @inline((l,n) -> σw*h/σf*(∂(r->∂(r->_s(l,n,r),r), r1) - l*(l+1)/r1^2*_s(l,n,r1)) + _s(l,n,r1)*l/r1 + ∂(r->_s(l,n,r),r1)*(1 + l*μr*h/r1)), 
          )
    return fs
end

#10.1103/PhysRevE.88.053010
@inline function bcs_t(b::Basis{ThinWall,Sphere}) 
    @inline _t = (l,n,r) -> r*t(Basis{ThinWall,Sphere}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (@inline((l,n) -> σw*h/σf*∂(r->_t(l,n,r),r1) + _t(l,n,r1)), )
    return fs
end


lpmax(b::Basis{ThinWall,Sphere}) = b.N
ltmax(b::Basis{ThinWall,Sphere}) = b.N

end