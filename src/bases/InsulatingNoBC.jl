module InsulatingBasisNoBC

using SparseArrays


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t
using ..Poly: ∂

export InsulatingNoBC

struct InsulatingNoBC end

InsulatingNoBC(N; kwargs...) = Basis{InsulatingNoBC}(; N, BC=NoBC(), kwargs...)

# @inline function t(::Type{Basis{InsulatingNoBC}}, l,m,n,r) 
#     fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n)))
#     return fac * r^l * (jacobi(n,0,l+1/2, 2r^2-1) - jacobi(n-1,0,l+1/2,2r^2-1)) 
# end

# @inline function s(::Type{Basis{InsulatingNoBC}}, l,m,n,r) 
#     fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
#     return fac * r^l * ( (2*l + 4*n - 3) * jacobi(n,0,l+1/2,2*r^2-1) - 2*(2*l + 4*n - 1)*jacobi(n-1,0,l+1/2,2*r^2-1) +(2*l + 4*n + 1)*jacobi(n-2,0,l+1/2,2*r^2-1))
# end

# @inline _nrange_p(b::Basis{InsulatingNoBC},l) = 1:((b.N-l+1)÷2)
# @inline _nrange_t(b::Basis{InsulatingNoBC},l) = 1:((b.N-l)÷2)

# @inline lpmax(b::Basis{InsulatingNoBC}) = b.N
# @inline ltmax(b::Basis{InsulatingNoBC}) = b.N


@inline function t(::Type{Basis{InsulatingNoBC}}, l, m, n, r)
    fac = sqrt(3 + 2l + 4n) / sqrt(l * (l + 1))
    return r^l * jacobi(n, 0, l + 1 / 2, 2r^2 - 1) * fac
end

@inline function s(::Type{Basis{InsulatingNoBC}}, l, m, n, r)
    fac = sqrt(3 + 2l + 4n) / sqrt(l * (l + 1))
    return r^l * jacobi(n, 0, l + 1 / 2, 2r^2 - 1) * fac
end

@inline _nrange_p(b::Basis{InsulatingNoBC}, l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{InsulatingNoBC}, l) = 0:((b.N-l)÷2)

@inline lpmax(b::Basis{InsulatingNoBC}) = b.N
@inline ltmax(b::Basis{InsulatingNoBC}) = b.N

@inline function bcs_t(::Type{Basis{InsulatingNoBC}})
    fs = (@inline((l, n) -> t(Basis{InsulatingNoBC}, l, 0, n, 1.0)),)
    return fs
end

@inline function bcs_p(::Type{Basis{InsulatingNoBC}})
    fs = (@inline((l, n) -> ∂(r -> s(Basis{InsulatingNoBC}, l, 0, n, r), 1.0) + (l + 1) * s(Basis{InsulatingNoBC}, l, 0, n, 1.0)),)
    return fs
end

end