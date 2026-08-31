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
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t, _lmn2cdeg_p, _lmn2cdeg_t

import ..Limace: inertial, _inertial_ss, _inertial_tt, _diffusion_ss, _diffusion_tt
import ..Quadrature: rquad
import ..Limace: inertial, diffusion

export ThinWall

struct ThinWall{q,h}; end

#based on appendix of Guervilly et al. (10.1103/PhysRevE.88.053010)

function ThinWall(N; q=1.0, h = 0.0, kwargs...)
    params=Dict(:q =>q, :h => h)
    return Basis{ThinWall{q,h},Sphere}(;N, V=Sphere(), BC=ThinWallBC{q,h}(), params,  kwargs...)
end

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

#based on Roberts et al. (2010)

# function ThinWall(N; c=0.0, c1 = 0.0, kwargs...)
#     params=Dict(:c =>c, :c1 => c1)
#     return Basis{ThinWall{c,c1},Sphere}(;N, V=Sphere(), BC=ThinWallBC{c,c1}(), params,  kwargs...)
# end

# function s(::Type{Basis{ThinWall{c,c1},Sphere}}, V::Sphere, l,m,n,r) where {c,c1}
#     fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n))) # ∫s⋅s dV ≠ 1 for h!=0.
#     coeff1 = (-3 + 2*l + 4*n)*(2 + 2*c1*l + c*(2 + c1*l)*(-1 + n)*(-1 + 2*l + 2*n)) 
#     coeff2 = -(-1 + 2*l + 4*n)*(4 + 4*c1*l + c*(2 + c1*l)*(3 - 2*n + 4*n^2 + l*(-2 + 4*n)))
#     coeff3 = (1 + 2*l + 4*n)*(2 + 2*c*n*(1 + 2*l + 2*n) + c1*l*(2 + c*n*(1 + 2*l + 2*n)))
#     return fac*r^l*(coeff1*jacobi(n,0,l+1/2,2r^2-1) + coeff2*jacobi(n-1,0,l+1/2,2r^2-1) + coeff3*jacobi(n-2,0,l+1/2,2r^2-1))
# end

# function t(::Type{Basis{ThinWall{c,c1},Sphere}}, V::Sphere, l,m,n,r) where {c,c1}
#     fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n))) # ∫t⋅t dV ≠ 1 for h!=0.
#     coeff1 = (1 + c*(l + n)*(-1 + 2*n))
#     coeff2 = -(1 + c*(1 + l + n)*(1 + 2*n))
#     return fac*r^l*(jacobi(n,0,l+1/2,2r^2-1) + coeff2/coeff1*jacobi(n-1,0,l+1/2,2r^2-1))
# end

@inline _nrange_p(b::Basis{ThinWall{q,h},Sphere},l) where {q,h} = 1:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{ThinWall{q,h},Sphere},l) where {q,h} = 1:((b.N-l)÷2)

lpmax(b::Basis{ThinWall{q,h},Sphere}) where {q,h} = b.N
ltmax(b::Basis{ThinWall{q,h},Sphere}) where {q,h} = b.N

_lmn2cdeg_p(b::Basis{ThinWall{q,h},Sphere}, l,m,n) where {q,h}  = l+2n-1
_lmn2cdeg_t(b::Basis{ThinWall{q,h},Sphere}, l,m,n) where {q,h} = l+2n


function diffusion(b::Basis{ThinWall{q,h}, Sphere}; η::T=1.0, threads=false, external=true) where {q,h,T}
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)
    r, wr = rquad(b.N+5,b.V)

    diff_s = [_diffusion_ss(b,b, (l,m,n), (l,m,n), r, wr)*η for (l,m,n) in lmnp]
    diff_t = [_diffusion_tt(b,b, (l,m,n), (l,m,n), r, wr)*η for (l,m,n) in lmnt]
    return spdiagm(vcat(diff_s,diff_t))
end

function inertial(b::Basis{ThinWall{q,h}, Sphere}; threads=false, external=true) where {q,h}
    T = typeof(b.V.r1)
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)

    r, wr = rquad(b.N+5,b.V)


    inertial_s = [T(_inertial_ss(b,b, (l,m,n), (l,m,n), r, wr)) for (l,m,n) in lmnp]
    inertial_s2 = [(n-1 > 0 ? T(_inertial_ss(b,b, (l,m,n), (l,m,n-1),r,wr)) : zero(T)) for (l,m,n) in lmnp[2:end]]
    inertial_t = [T(_inertial_tt(b,b, (l,m,n), (l,m,n), r, wr)) for (l,m,n) in lmnt]
    inertial_t2 = [(n-1 > 0 ? T(_inertial_tt(b,b, (l,m,n), (l,m,n-1),r,wr)) : zero(T)) for (l,m,n) in lmnt[2:end]]

    d = vcat(inertial_s,inertial_t)
    d2 = [inertial_s2;0.0; inertial_t2]
    return complex(sparse(SymTridiagonal(d, d2)))
end 





end