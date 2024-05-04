module ViscousBasis


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s
import ..Limace: _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!
import ..Limace: inertial, diffusion

export Viscous

struct Viscous; end

Viscous(N; kwargs...) = Basis{Viscous}(;N, BC=NoSlipBC(), V=Sphere(), kwargs...)

"""
$(TYPEDSIGNATURES)

Chen et al. (2018) toroidal scalar, orthogonal w.r.t ∫ u⋅∇²u dV with 0 ≤ r ≤ 1.
"""
@inline function t(::Type{Basis{Viscous}}, V::Volume, l,m,n,r) 
    fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n)))
    return fac * r^l * (jacobi(n,0,l+1/2, 2r^2-1) - jacobi(n-1,0,l+1/2,2r^2-1)) 
end

"""
$(TYPEDSIGNATURES)

Chen et al. (2018) (2.38), (2.39) poloidal scalar, orthogonal w.r.t ∫ u⋅∇²u dV with 0 ≤ r ≤ 1.
"""
@inline function s(::Type{Basis{Viscous}}, V::Volume, l,m,n,r) 
    c1 = 2l+4n+1
    c2 = -2(2l+4n+3)
    c3 = 2l+4n+5
    fac = 1/(sqrt(2l*(1 + l)*(1 + 2*l + 4*n)*(3 + 2*l + 4*n)*(5 + 2*l + 4*n)))
    return fac*r^l*(c1*jacobi(n+1,0,l+1/2,2r^2-1) + c2*jacobi(n,0,l+1/2,2r^2-1) + c3*jacobi(n-1,0,l+1/2,2r^2-1)  ) 
end

@inline _nrange_p(b::Basis{Viscous},l) = 1:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{Viscous},l) = 1:((b.N-l)÷2+1)

@inline lpmax(b::Basis{Viscous}) = b.N-1
@inline ltmax(b::Basis{Viscous}) = b.N


#Inertial term/inner products

@inline function _inertial_tt(l,n,n2)
    if n==n2 
        return one(l)
    elseif (n==n2+1)
        return -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n)) + 3/(5 + 2*l + 4*(-1 + n)))/2
    elseif (n==n2-1)
        return -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n2)) + 3/(5 + 2*l + 4*(-1 + n2)))/2
    else
        return zero(l)
    end
end

@inline function _inertial_ss(l,n,n2)
    if n==n2 
        return one(l)
    elseif (n==n2+1) 
        return -1/2sqrt(1 + 3/(2 + 4*l + 8*(-1 + n)) - 3/(18 + 4*l + 8*(-1 + n)))
    elseif (n==n2-1)
        return -1/2sqrt(1 + 3/(2 + 4*l + 8*(-1 + n2)) - 3/(18 + 4*l + 8*(-1 + n2)))
    else 
        return zero(l)
    end
end

function inertial(b::Basis{Viscous}, ::Type{T}=Float64) where {T<:Number}
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)


    inertial_s = [_inertial_ss(T(l), T(n), T(n)) for (l,m,n) in lmnp]
    inertial_s2 = [(n-1 > 0 ? _inertial_ss(T(l), T(n), T(n-1)) : 0.0) for (l,m,n) in lmnp[2:end]]
    inertial_t = [_inertial_tt(T(l), T(n), T(n)) for (l,m,n) in lmnt]
    inertial_t2 = [(n-1 > 0 ? _inertial_tt(T(l), T(n), T(n-1)) : 0.0) for (l,m,n) in lmnt[2:end]]

    d = vcat(inertial_s,inertial_t)
    d2 = [inertial_s2;0.0; inertial_t2]
    return SymTridiagonal(d, d2)
end 

#Coriolis force

@inline function _coriolis_tt(l,m,n,n2; Ω::T = 2.0) where T
    if n==n2
        return Ω*im*m/(l*(l+1)) #factor 2??
    elseif (n==n2+1)
        n = n-1
        return -Ω*im*m/(2l*(1 + l)*sqrt(1 + 3/(2*(-1 + 2*l + 4*n)) - 3/(14 + 4*l + 8*n)))
    elseif (n==n2-1)
        n = n+1
        return -Ω*im*m*sqrt(1 - 3/(-3 + 2*l + 4*n) + 3/(1 + 2*l + 4*n))/(2l*(1 + l))
    end
    return zero(T)
end

@inline function _coriolis_ss(l,m,n,n2; Ω::T = 2.0) where T 
    if n==n2
        return Ω*im*m/(l*(l+1)) #factor 2??
    elseif (n==n2+1)
        n = n-1
        return -Ω*im*m/(2l*(1 + l)*sqrt(1 + 3/(2 + 4*l + 8*n) - 3/(18 + 4*l + 8*n)))
    elseif (n==n2-1)
        n=n+1
        return -Ω*im*m*sqrt(1 - 3/(-1 + 2*l + 4*n) + 3/(3 + 2*l + 4*n))/(2l*(1 + l))    
    end
    return zero(T)
end


@inline function _coriolis_ts(l,l2,m,m2,n,n2; Ω::T = 2.0) where T
    @assert m==m2
    if (l==l2+1)
        if n==n2
            return Ω*((-1 + l^2)*sqrt(((l - m)*(l + m))/(1 - 5*l^2 + 4*l^4)))/l
        elseif (n == n2-1)
            return -Ω*((-1 + l^2)*sqrt(((l - m)*(l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((1 - 5*l^2 + 4*l^4)*(1 + 2*l + 4*n)*(5 + 2*l + 4*n))))/l/2 
        elseif (n == n2+1)
            return -Ω*((-1 + l^2)*sqrt((l - m)*(l + m))*sqrt(-5 + 2*l + 4*n))/(l*sqrt(((1 - 5*l^2 + 4*l^4)*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
        end
    elseif (l==l2-1)
        if n==n2
            return -Ω*sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((3 + 4*l*(2 + l))*(1 + 2*l + 4*n)*(5 + 2*l + 4*n)))/(1 + l)/2
        elseif n==n2+1
            return Ω*sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m))/(3 + 4*l*(2 + l)))/(1 + l)
        elseif n==n2+2
            return -Ω*(sqrt(l*(2 + l)*(1 + l - m)*(1 + l + m))*sqrt(-5 + 2*l + 4*n))/((1 + l)*sqrt(((3 + 4*l*(2 + l))*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
        else
        end
    end
    return zero(T)
end

function _coriolis_st(l,l2,m,m2,n,n2; Ω::T = 2.0) where T
    aij = -_coriolis_ts(l2,l,m2,m,n2,n; Ω)
    return aij
end
# @inline function _coriolis_st(l2,l,m2,m,n2,n; Ω::T = 2.0) where T
#     @assert m==m2
#     if (l==l2+1)
#         if n==n2
#             return -Ω*((-1 + l^2)*sqrt(((l - m)*(l + m))/(1 - 5*l^2 + 4*l^4)))/l
#         elseif (n == n2-1)
#             return Ω*((-1 + l^2)*sqrt(((l - m)*(l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((1 - 5*l^2 + 4*l^4)*(1 + 2*l + 4*n)*(5 + 2*l + 4*n))))/l/2 
#         elseif (n == n2+1)
#             return Ω*((-1 + l^2)*sqrt((l - m)*(l + m))*sqrt(-5 + 2*l + 4*n))/(l*sqrt(((1 - 5*l^2 + 4*l^4)*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
#         end
#     elseif (l==l2-1)
#         if n==n2
#             return Ω*sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((3 + 4*l*(2 + l))*(1 + 2*l + 4*n)*(5 + 2*l + 4*n)))/(1 + l)/2
#         elseif n==n2+1
#             return -Ω*sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m))/(3 + 4*l*(2 + l)))/(1 + l)
#         elseif n==n2+2
#             return Ω*(sqrt(l*(2 + l)*(1 + l - m)*(1 + l + m))*sqrt(-5 + 2*l + 4*n))/((1 + l)*sqrt(((3 + 4*l*(2 + l))*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
#         end
#     else
#         return zero(T)
#     end
# end

function _coriolis_poloidal_poloidal!(b::Basis{Viscous}, is, js, aijs, lmn2k_p, l, m, r, wr, Ω::T; applyBC=true) where T
    for n in nrange_p(b,l)
        njs_all = nrange_p(b,l)
        for n2 in max(n-1,first(njs_all)):min(n+1, last(njs_all))
            aij = _coriolis_ss(T(l), T(m), T(n), T(n2); Ω)
            appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
        end
    end
    return nothing
end

function _coriolis_poloidal_toroidal!(b::Basis{Viscous}, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω::T) where T
    for n in nrange_p(b,l)
        njs_all = nrange_t(b,l2)
        for n2 in max(n-2,first(njs_all)):min(n+2, last(njs_all))
            aij = _coriolis_st(T(l),T(l2),T(m),T(m),T(n),T(n2); Ω)
            appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_t[(l2,m,n2)] + _np, aij)
        end
    end
    return nothing
end


function _coriolis_toroidal_toroidal!(b::Basis{Viscous}, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω::T; applyBC=true) where T
    for n in nrange_t(b,l)
        njs_all = nrange_t(b,l)
        for n2 in max(n-1,first(njs_all)):min(n+1, last(njs_all))
            aij = _coriolis_tt(T(l), T(m), T(n), T(n2); Ω)
            appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
        end
    end
    return nothing
end

function _coriolis_toroidal_poloidal!(b::Basis{Viscous}, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω::T) where T
    for n in nrange_t(b,l) 
        njs_all = nrange_p(b,l2)
        for n2 in max(n-2,first(njs_all)):min(n+2, last(njs_all))
            aij = _coriolis_ts(T(l),T(l2),T(m),T(m),T(n),T(n2); Ω)
            appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_p[(l2,m,n2)], aij)
        end
    end
    return nothing
end


#Viscosity


@inline function _diffusion_tt(l,n; ν::T = 1.0) where T
    return -ν*((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
end


@inline function _diffusion_ss(l,n; ν::T = 1.0) where T
    return -ν*((1 + 2*l + 4*n)*(5 + 2*l + 4*n))/2
end

function diffusion(b::Basis{Viscous}; ν::T=1.0, applyBC=true) where T
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)

    diff_s = [_diffusion_ss(T(l), T(n); ν) for (l,m,n) in lmnp]
    diff_t = [_diffusion_tt(T(l), T(n);ν) for (l,m,n) in lmnt]
    return spdiagm(vcat(diff_s,diff_t))
end


end