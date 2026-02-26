module PerfectlyConductingBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, _lmn2cdeg_p, _lmn2cdeg_t
import ..Limace: inertial, _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!

export PerfectlyConducting 

struct PerfectlyConducting; end

PerfectlyConducting(N; kwargs...) = Basis{PerfectlyConducting, Sphere}(;N, BC=PerfectlyConductingBC(), V=Sphere(), kwargs...)

"""
    t(::Type{Basis{PerfectlyConducting, Sphere}}, V::Sphere, l,m,n,r)

```math
t_{l,n,m}(r) = f_{l,n}r^l J_n^{(0,l+1/2)}(2r^2-1)
```
with 

```math
f_{l,n} = \\sqrt{\\frac{3+2l+4n}{l(l+1)}}
```

[livermore_compendium_2014](@citet) (5.1), normalized to unit energy.
"""
@inline function t(::Type{Basis{PerfectlyConducting, Sphere}}, V::Sphere, l,m,n,r)
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

"""
    s(::Type{Basis{PerfectlyConducting, Sphere}}, V::Sphere, l,m,n,r)

```math
s_{l,n,m}(r) = f_{l,n}(1-r^2)r^l J_n^{(1,l+1/2)}(2r^2-1)
```

with
```math
f_{l,n} = \\sqrt{\\frac{5+2l+4n}{4l(l+1)(n+1)^2}}
```

[livermore_compendium_2014](@citet) (5.6), normalized to unit energy. 
"""
@inline function s(::Type{Basis{PerfectlyConducting, Sphere}}, V::Sphere, l,m,n,r)
    fac = sqrt(5+2l+4n)/sqrt(4l*(l+1)*(n+1)^2)
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end

@inline _nrange_p(b::Basis{PerfectlyConducting, Sphere},l) = 0:((b.N-l+1)÷2-1)
@inline _nrange_t(b::Basis{PerfectlyConducting, Sphere},l) = 0:((b.N-l)÷2)

@inline lpmax(b::Basis{PerfectlyConducting, Sphere}) = b.N
@inline ltmax(b::Basis{PerfectlyConducting, Sphere}) = b.N

n(N) = (2N^3+9N^2+7N)÷6

_np(N::Int) = ((-1)^N*(3 + (-1)^N*(-3+2N*(-1+N*(3+N)))))÷12

@inline function _np(N::Int, m::Int) 
    if (m == 0)
        m+=1
    end
    return (N-abs(m)+1)^2 ÷ 4
end

@inline function np(b::Basis{PerfectlyConducting, Sphere})
    if isaxisymmetric(b)
        return _np(b.N,first(b.m))
    else
        return _np(b.N)
    end
end

_nt(N) = ((-1)^N*(-3 + (-1)^N*(3 + 2*N*(2 + N)*(4 + N))))÷12

@inline function _nt(N::Int, m::Int)
    if (m==0)
        m+=1
    end
    return (N-abs(m)+2)^2 ÷ 4
end

@inline function nt(b::Basis{PerfectlyConducting, Sphere})
    if isaxisymmetric(b)
        return _nt(b.N,first(b.m))
    else
        return _nt(b.N)
    end
end


function inertial(b::Basis{PerfectlyConducting, Sphere}; kwargs...)
    return one(typeof(b.V.r1))*I(length(b))
end

end