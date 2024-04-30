module UnconstrainedBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t

export Unconstrained

struct Unconstrained; end

Unconstrained(N; kwargs...) = Basis{Unconstrained}(;N, BC=NoBC(), kwargs...)

@inline function t(::Type{Basis{Unconstrained}}, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return fac*r^l*jacobi(n,0,l+1/2, 2r^2-1)
end

@inline function s(::Type{Basis{Unconstrained}}, l,m,n,r)
    fac = n==1 ? 1/sqrt(l*(l+1)^2) : 1/sqrt(l*(l+1)*(-3+2l+4n))
    return fac*r^l*(jacobi(n-1,0,l+1/2, 2r^2-1) - jacobi(n-2,0,l+1/2, 2r^2-1))
end

@inline _nrange_p(b::Basis{Unconstrained},l) = 1:((b.N-l+1)÷2+1)
@inline _nrange_t(b::Basis{Unconstrained},l) = 0:((b.N-l)÷2)

# @inline nrange_p_bc(b::Basis{Unconstrained},l) = 0:((b.N-l+1)÷2-1)
# @inline nrange_t_bc(b::Basis{Unconstrained},l) = nrange_t(b, l)

@inline function bcs_p(::Type{Basis{Unconstrained}})
    return ()
end

@inline function bcs_t(::Type{Basis{Unconstrained}})
    return ()
end


lpmax(b::Basis{Unconstrained}) = b.N-1
ltmax(b::Basis{Unconstrained}) = b.N


end