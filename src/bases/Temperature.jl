module TemperatureBasis

using SparseArrays


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, _lmn2cdeg_p, _lmn2cdeg_t
import ..Limace: inertial, diffusion

export Temperature

struct Temperature; end

Temperature(N; kwargs...) = Basis{Temperature, Sphere}(;N, BC=DirichletBC(), V=Sphere(), kwargs...)

"""
    t(::Type{Basis{Temperature, Sphere}}, V::Sphere, l,m,n,r)

```math
t_{l,m,n}(r) = f_{l,n}r^l\\left( J_{n}^{(0,l+1/2)}(2r^2-1) - J_{n-1}^{(0,l+1/2)}(2r^2-1)\\right)
```
with
```math
f_{l,n} = \\left(1/(2l+4n-1) + 1/(2l+4n+3) \\right)^{-1/2}
```

[gerick_interannual_2024](@citet) (A7)
"""
@inline function t(::Type{Basis{Temperature, Sphere}}, V::Sphere, l,m,n,r) 
    fac = 1/sqrt(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n))
    return fac * r^l * (jacobi(n,0,l+1/2, 2r^2-1) - jacobi(n-1,0,l+1/2,2r^2-1)) 
end

@inline _nrange_p(b::Basis{Temperature, Sphere},l) = 0:-1
@inline _nrange_t(b::Basis{Temperature, Sphere},l) = 1:((b.N-l)÷2)

@inline lpmax(b::Basis{Temperature, Sphere}) = 0
@inline ltmax(b::Basis{Temperature, Sphere}) = b.N


_np(N) = 0
_nt(N) = ((-1)^N*(-3 + (-1)^N*(3 - 8*N + 2*N^3)))÷12 #+ 1

_nlp(N,l) = 0
_nlt(N,l) = ((-1)^N*(-3 + 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(-11 - 4*l*(3 + l) + 12*N + 6*l*N))))÷12

lmn2k_p(l,m,n,N) = 0
lmn2k_t(l,m,n,N) = _nlt(N,l-1) + (l+m)*((N-l)÷2) + n

_lmn2cdeg_p(b::Basis{Temperature, Sphere}, l,m,n) = 0
_lmn2cdeg_t(b::Basis{Temperature, Sphere}, l,m,n) = l+2n


# following equaitons have a 1/(l*(l+1)) factor in inertial() and diffusion(), to compensate for ∫ TᵢTⱼ dV = 1/(l(l+1)).
# This is because 
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

function inertial(b::Basis{Temperature, Sphere}; threads=false, external=true)
    T = typeof(b.V.r1)
    lmnt = lmn_t(b)


    inertial_t = [_inertial_tt(T(l), T(n), T(n)) for (l,m,n) in lmnt]
    inertial_t2 = [(n-1 > 0 ? _inertial_tt(T(l), T(n), T(n-1)) : 0.0) for (l,m,n) in lmnt[2:end]]

    return SymTridiagonal(inertial_t, inertial_t2)
end 

#diffusion

@inline function _diffusion_tt(l,n; η::T = 1.0) where T 
    return -η*((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
end

function diffusion(b::Basis{Temperature, Sphere}; η::T=1.0, threads=false, external=true) where T
    lmnt = lmn_t(b)


    diff_t = [_diffusion_tt(T(l), T(n); η)  for (l,m,n) in lmnt]
    return spdiagm(diff_t)
end

end