module ThinWallBasis

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

import ..Limace: inertial, _inertial_ss, _inertial_tt
import ..Quadrature: rquad

export ThinWall

struct ThinWall; end

function ThinWall(N; σw=1.0, σf = 1.0, h = 0.0, kwargs...)
    params=Dict(:σw => σw, :σf => σf, :h => h)
    return Basis{ThinWall}(;N, V=Sphere(), BC=NoBC(), params,  kwargs...)
end

s(::Type{Basis{ThinWall}}, V::Volume, l,m,n,r) = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{ThinWall}}, V::Volume, l,m,n,r) = t(Basis{Unconstrained}, V, l,m,n,r) 

@inline _nrange_p(b::Basis{ThinWall},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{ThinWall},l) = 0:((b.N-l)÷2)

#10.1103/PhysRevE.88.053010
@inline function bcs_p(b::Basis{ThinWall}) 
    @inline _s = (l,n,r) -> r*s(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (
          @inline((l,n) -> σw*h/σf*(∂(r->∂(r->_s(l,n,r),r), r1) - l*(l+1)/r1^2*_s(l,n,r1)) + _s(l,n,r1)*l/r1 + ∂(r->_s(l,n,r),r1)*(1 + l*h/r1)), 
          )
    return fs
end

#10.1103/PhysRevE.88.053010
@inline function bcs_t(b::Basis{ThinWall}) 
    @inline _t = (l,n,r) -> t(Basis{ThinWall}, b.V, l, 0, n, r)
    (; r1) = b.V 
    h, σf, σw = b.params[:h], b.params[:σf], b.params[:σw]
    fs = (@inline((l,n) -> σw/σf*h*∂(r->_t(l,n,r),r1) + _t(l,n,r1)), )
    return fs
end


lpmax(b::Basis{ThinWall}) = b.N
ltmax(b::Basis{ThinWall}) = b.N

function inertial(b::Basis{ThinWall}, ::Type{T}=Float64; external=false) where {T<:Number}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    nu = length(b)

    #m == m2 and only l==l2 needs to be considered.
    @inbounds for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            for n in nrange_p_bc(b,l)
                aij = _inertial_ss(b, (l,m,n), (l,m,n), r,wr; external)
                appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n)], aij)
            end
        end
    end

    @inbounds for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            for n in nrange_t_bc(b,l)
                aij = _inertial_tt(b, (l,m,n), (l,m,n), r,wr)
                appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n)] + _np, aij)
            end
        end
    end


    return sparse(is, js, aijs, nu, nu)
end


end