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
import ..Limace: inertial, diffusion

export ThinWall

struct ThinWall{q,h}; end

function ThinWall(N; q=1.0, h = 0.0, kwargs...)
    params=Dict(:q =>q, :h => h)
    return Basis{ThinWall{q,h},Sphere}(;N, V=Sphere(), BC=ThinWallBC{q,h}(), params,  kwargs...)
end

#based on appendix of 10.1103/PhysRevE.88.053010
function s(::Type{Basis{ThinWall{q,h},Sphere}}, V::Sphere, l,m,n,r) where {q,h}
    fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n))) # ∫s⋅s dV ≠ 1 for h!=0.
    coeff1 = (-3 + 2*l + 4*n)*(1 + h*(l + 2*l*(-1 + n)*q + (1 - 3*n + 2*n^2)*q))
    coeff2 = -(-1 + 2*l + 4*n)*(2 + h*(3 - 2*n + 4*n^2)*q + h*l*(2 + (-2 + 4*n)*q))
    coeff3 = (1 + 2*l + 4*n)*(1 + h*(l + 2*l*n*q + n*(1 + 2*n)*q)) 
    return fac*r^l*(coeff1*jacobi(n,0,l+1/2,2r^2-1) + coeff2*jacobi(n-1,0,l+1/2,2r^2-1) + coeff3*jacobi(n-2,0,l+1/2,2r^2-1))
end

function t(::Type{Basis{ThinWall{q,h},Sphere}}, V::Sphere, l,m,n,r) where {q,h}
    fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n))) # ∫t⋅t dV ≠ 1 for h!=0.
    coeff1 = (1 + h*(l + n)*(-1 + 2*n)*q)
    coeff2 = -(1 + h*(1 + l + n)*(1 + 2*n)*q)
    return fac*r^l*(jacobi(n,0,l+1/2,2r^2-1) + coeff2/coeff1*jacobi(n-1,0,l+1/2,2r^2-1))
end

@inline _nrange_p(b::Basis{ThinWall{q,h},Sphere},l) where {q,h} = 1:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{ThinWall{q,h},Sphere},l) where {q,h} = 1:((b.N-l)÷2)

lpmax(b::Basis{ThinWall{q,h},Sphere}) where {q,h} = b.N
ltmax(b::Basis{ThinWall{q,h},Sphere}) where {q,h} = b.N


end