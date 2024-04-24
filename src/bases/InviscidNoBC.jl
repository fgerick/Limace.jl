module InviscidBasisNoBC

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t
import ..Limace: _coriolis_poloidal_poloidal!, _coriolis_toroidal_toroidal!,  _coriolis_poloidal_toroidal!, _coriolis_toroidal_poloidal!

export InviscidNoBC

struct InviscidNoBC; end

InviscidNoBC(N; kwargs...) = Basis{InviscidNoBC}(;N, BC=NoBC(), kwargs...)

@inline function t(::Type{Basis{InviscidNoBC}}, l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end

@inline function s(::Type{Basis{InviscidNoBC}}, l,m,n,r)
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
    # return r^l*ultrasphericalc(n,l-1//2,2r^2-1)
    # return r^l*(jacobi(n+1,0,l+1/2, 2r^2-1) - jacobi(n,0,l+1/2, 2r^2-1))
end

@inline _nrange_p(b::Basis{InviscidNoBC},l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{InviscidNoBC},l) = 0:((b.N-l)÷2)

@inline nrange_p_bc(b::Basis{InviscidNoBC},l) = 0:((b.N-l+1)÷2-1)
@inline nrange_t_bc(b::Basis{InviscidNoBC},l) = nrange_t(b, l)

@inline function bcs_p(::Type{Basis{InviscidNoBC}})
    fs = (@inline((l,n)->s(Basis{InviscidNoBC}, l, 0, n, 1.0)), )
    return fs
end

@inline function bcs_t(::Type{Basis{InviscidNoBC}})
    return ()
end


lpmax(b::Basis{InviscidNoBC}) = b.N-1
ltmax(b::Basis{InviscidNoBC}) = b.N


end