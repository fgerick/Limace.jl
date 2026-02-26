module ViscousShellBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t, SphericalShell
import ..ViscousBasis: Viscous
export ViscousShell

ViscousShell(N; r0 = 0.35, r1=1.0, kwargs...) = Basis{Viscous, SphericalShell}(;N, BC=NoBC(), V=SphericalShell(r0,r1), kwargs...)

@inline function t(::Type{Basis{Viscous, SphericalShell}}, V::Volume, l::Int,m::Int,n::Int,r::T) where T 
	r0,r1 = V.r0, V.r1
	x = (2r-(r1+r0))/(r1-r0) #map to x ∈ [-1,1]
	return (r0-r)*(r1-r)*jacobi(n-2,0.0,0.0, x) # - jacobi(n-1,0.0,0.0, x)
end

@inline function s(::Type{Basis{Viscous, SphericalShell}}, V::Volume, l::Int,m::Int,n::Int,r::T) where T
	r0,r1 = V.r0, V.r1
	x = (2r-(r1+r0))/(r1-r0) #map to x ∈ [-1,1]
	return (r0-r)*(r1-r)*jacobi(n-2,0.0,0.0, x) #- jacobi(n-1,0.0,0.0, x)
end

@inline _nrange_p(b::Basis{Viscous, SphericalShell},l) = 2:((b.N-l+1)+3)
@inline _nrange_t(b::Basis{Viscous, SphericalShell},l) = 2:((b.N-l)+2)

# @inline nrange_p_bc(b::Basis{ViscousShell},l) = 0:((b.N-l+1)÷2-1)
# @inline nrange_t_bc(b::Basis{ViscousShell},l) = nrange_t(b, l)

@inline function bcs_p(b::Basis{Viscous, SphericalShell})
    fs = (@inline((l, n) -> ∂(r -> s(Basis{Viscous, SphericalShell}, b.V, l, 0, n, r), b.V.r0)), 
		  @inline((l, n) -> ∂(r -> s(Basis{Viscous, SphericalShell}, b.V, l, 0, n, r), b.V.r1)),)
    return fs
end

@inline function bcs_t(b::Basis{Viscous, SphericalShell})
    fs = ()
	return fs
end


lpmax(b::Basis{Viscous, SphericalShell}) = b.N
ltmax(b::Basis{Viscous, SphericalShell}) = b.N


end