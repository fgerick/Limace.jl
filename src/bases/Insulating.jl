module InsulatingBasis

using SparseArrays


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s
import ..Limace: inertial, diffusion

export Insulating

struct Insulating; end

Insulating(N; kwargs...) = Basis{Insulating}(;N, BC=InsulatingBC(), V=Sphere(), kwargs...)

@inline function t(::Type{Basis{Insulating}}, V::Volume, l,m,n,r) 
    fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n)))
    return fac * r^l * (jacobi(n,0,l+1/2, 2r^2-1) - jacobi(n-1,0,l+1/2,2r^2-1)) 
end

@inline function s(::Type{Basis{Insulating}}, V::Volume, l,m,n,r) 
    fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    return fac * r^l * ( (2*l + 4*n - 3) * jacobi(n,0,l+1/2,2*r^2-1) - 2*(2*l + 4*n - 1)*jacobi(n-1,0,l+1/2,2*r^2-1) +(2*l + 4*n + 1)*jacobi(n-2,0,l+1/2,2*r^2-1))
end

@inline _nrange_p(b::Basis{Insulating},l) = 1:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{Insulating},l) = 1:((b.N-l)÷2)

@inline lpmax(b::Basis{Insulating}) = b.N
@inline ltmax(b::Basis{Insulating}) = b.N


n(N) = ((-1)^(2*N)*(-1 + N)*N*(5 + 2*N))÷6
_np(N) = ((-1)^N*(3 + (-1)^N*(-3 + 2*N*(-1 + N*(3 + N)))))÷12
_nt(N) = ((-1)^N*(-3 + (-1)^N*(3 - 8*N + 2*N^3)))÷12

_nlp(N,l) = ((-1)^N*(3 - 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(1 - 6*l - 4*l^2 + 6*(2 + l)*N))))÷12
_nlt(N,l) = ((-1)^N*(-3 + 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(-11 - 4*l*(3 + l) + 12*N + 6*l*N))))÷12

lmn2k_p(l,m,n,N) = _nlp(N,l-1) + (l+m)*((N-l+1)÷2) + n
lmn2k_t(l,m,n,N) = _nlt(N,l-1) + (l+m)*((N-l)÷2) + n




#inner products

@inline function _inertial_tt(l,n,n2)
    if n==n2 
        return one(l)
    elseif (n==n2+1)
        return -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n)) + 3/(5 + 2*l + 4*(-1 + n)))/2
    elseif (n==n2-1)
        return -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n2)) + 3/(5 + 2*l + 4*(-1 + n2)))/2
    end
    return zero(l)
end

@inline function _inertial_ss(l,n,n2)
    if n==n2 
        return one(l)
    elseif (n==n2+1) 
        return -sqrt(1 + 3/(5 - 2*l - 4*n) + 3/(-1 + 2*l + 4*n))/2
    elseif (n==n2-1)
        return -sqrt(1 + 3/(5 - 2*l - 4*n2) + 3/(-1 + 2*l + 4*n2))/2
    end
    return zero(l)
end

function inertial(b::Basis{Insulating}, ::Type{T}=Float64; external=true) where {T<:Number}
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

# function inertial(b::Basis{Insulating}, ::Type{T}=Float64) where {T<:Number}
    
#     is, js, aijs = Int[], Int[], Complex{T}[]
#     lmn2k_p = lmn2k_p_dict(b)
#     lmn2k_t = lmn2k_t_dict(b)
#     _np = np(b)
#     nmat = length(b)


#     @inbounds for l in 1:lpmax(b)
#         for m in intersect(b.m, -l:l)
#             ns = nrange_p(b,l)
#             for n in ns
#                 for n2 in max(n-1,1):min(n+1,last(ns)) #symtridiagonal
#                     aij = _inertial_ss(T(l), T(n), T(n2))
#                     appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
#                 end
#             end
#         end
#     end

#     @inbounds for l in 1:ltmax(b)
#         for m in intersect(b.m, -l:l)
#             ns = nrange_t(b,l)
#             for n in ns
#                 for n2 in max(n-1,1):min(n+1,last(ns)) #symtridiagonal
#                     aij = _inertial_tt(T(l), T(n), T(n2))
#                     appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
#                 end
#             end
#         end
#     end

#     return SymTridiagonal(sparse(is,js,aijs, nmat, nmat))
# end

#diffusion

@inline function _diffusion_tt(l,n; η::T = 1.0) where T 
    return -η*((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
end

@inline function _diffusion_ss(l,n; η::T = 1.0) where T
    return -η*((-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/2
end

function diffusion(b::Basis{Insulating}; η::T=1.0, applyBC=true, external=true) where T
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)


    diff_s = [_diffusion_ss(T(l), T(n); η) for (l,m,n) in lmnp]
    diff_t = [_diffusion_tt(T(l), T(n); η) for (l,m,n) in lmnt]
    return spdiagm(vcat(diff_s,diff_t))
end


"""
$(TYPEDSIGNATURES)

The basis elements are normalized so that ∫₀^∞ B⋅B dV = 1.
There is no external contribution for the toroidal components. 
For the poloidal components with n=1, there is an external contribution (r>1). 
To normalize the basis elements so that ∫₀¹ B⋅B dV = 1, the following norm function can be used:
"""
function unitspherenorm(l,n)
    if n==1
        return sqrt((6+8l*(l+2))/(6+l*(6l+11)))
    else
        return 1.0
    end
end

"""
$(TYPEDSIGNATURES)

Convenience functions to normalize a B₀ that is assembled from different components to have (4π/3)⁻¹ ∫₀¹ B⋅B dV = 1
Construct Aᵢⱼ = ∫₀¹ Bᵢ ⋅Bⱼ dV 
"""
function inner_b0norm(lmn_p, lmn_t)
    T= Float64

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        _inner_ss!(is,js,aijs,i,i, T(l),T(n),T(n))
        aijs[end]/=unitspherenorm(l,n)^2
        if n>1
            _inner_ss!(is,js,aijs,i,i-1, T(l),T(n),T(n-1))
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        _inner_tt!(is,js,aijs, i+np,i+np, T(l),T(n),T(n))
        if n>1
            _inner_tt!(is,js,aijs, i+np,i-1+np, T(l),T(n),T(n-1))
        end
    end

    LHS = sparse(is,js,aijs, nu, nu)
    return LHS

end

"""
$(TYPEDSIGNATURES)

Normalize `B0fac` (vector of scalars) corresponding to `lmn_p` and `lmn_t` basis elements. `length(B0fac)=length(lmn_p)+length(lmn_t)`.
Final array of `B0fac` will ensure that (4π/3)⁻¹ ∫₀¹ B⋅B dV = 1.
"""
function norm_B0fac!(B0fac, lmn_p, lmn_t)
    A = inner_b0norm(lmn_p, lmn_t)
    for (i,(l,m,n)) in enumerate(lmn_p)
        B0fac[i] *= unitspherenorm(l,n)
    end
    _n = sqrt(B0fac'*A*B0fac)
    B0fac.*=√(4π/3)/_n
    return nothing
end	




end