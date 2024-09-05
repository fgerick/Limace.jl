module UnconstrainedBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export Unconstrained

struct Unconstrained; end

Unconstrained(N; r=1.0, kwargs...) = Basis{Unconstrained}(;N, BC=NoBC(), V=Sphere(; r), kwargs...)

@inline function t(::Type{Basis{Unconstrained}}, V::Volume, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    x = (2r^2-V.r1^2)/V.r1^2
    return fac*r^l*jacobi(n,0,l+1/2, x)
end

@inline function s(::Type{Basis{Unconstrained}}, V::Volume, l,m,n,r)
    fac = 1/sqrt(l*(l+1)*(1+2l+4n))
    x = (2r^2-V.r1^2)/V.r1^2
    return fac*r^l*(jacobi(n,0,l+1/2, x) - jacobi(n-1,0,l+1/2, x))
end

@inline _nrange_p(b::Basis{Unconstrained},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{Unconstrained},l) = 0:((b.N-l)÷2)


@inline function bcs_p(b::Basis{Unconstrained})
    return ()
end

@inline function bcs_t(b::Basis{Unconstrained})
    return ()
end


lpmax(b::Basis{Unconstrained}) = b.N
ltmax(b::Basis{Unconstrained}) = b.N


end