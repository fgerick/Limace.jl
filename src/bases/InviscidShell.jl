module InviscidShellBasis

using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t, SphericalShell
import ..InviscidBasis: Inviscid
export InviscidShell


InviscidShell(N; r0 = 0.35, r1=1.0, kwargs...) = Basis{Inviscid, SphericalShell}(;N, BC=NoBC(), V=SphericalShell(r0,r1), kwargs...)

@inline function t(::Type{Basis{Inviscid, SphericalShell}}, V::Volume, l,m,n,r) 
	r0,r1 = V.r0, V.r1
	x = (2r-(r1+r0))/(r1-r0) #map to x ∈ [-1,1]
	fac = 1/sqrt(-l*(l+1)*(r0-r1)*((-2+3n+3n^2)*r0^2+2*(-1+n+n^2)*r0*r1+(-2+3n+3n^2)*r1^2)/(-6-4n+24n^2+16n^3))
	return fac*jacobi(n,0,0, x) # - jacobi(n-1,0.0,0.0, x)
end

# @inline function s(::Type{Basis{Inviscid, SphericalShell}}, V::Volume, l,m,n,r)
# 	r0,r1 = V.r0, V.r1
# 	x = (2r-(r1+r0))/(r1-r0) #map to x ∈ [-1,1]
# 	if n == 1
# 		return 1-x
# 	else
# 		α = -1/2
# 		β = 1.0
# 		c1	= n*(β + α + 2n)*(α + β + n + 2)		
# 		c2	= -(β + 1 + 2n + α)*(2n^2 + 2n*β + 2n + 2n*α + β*α + α^2 + 3α + 2 + β)		
# 		c3	= (α + n + 1)*(β + n - 1)*(β + α + 2n + 2)
# 		# return jacobi(n,-1/2,1.0, x) #- jacobi(n-1,3/2,1.0, x)
# 		return c1*jacobi(n,α,β,x) + c2*jacobi(n-1,α,β,x) + c3*jacobi(n-2,α,β,x)
# 	end
# end
# @inline _nrange_p(b::Basis{Inviscid, SphericalShell},l) = 1:(b.N-l)+1

@inline function s(::Type{Basis{Inviscid, SphericalShell}}, V::Volume, l,m,n,r) 
	r0,r1 = V.r0, V.r1
	x = (2r-(r1+r0))/(r1-r0) #map to x ∈ [-1,1]
	# fac = 1/sqrt(-(r0-r1)*((-2+3n+3n^2)*r0^2+2*(-1+n+n^2)*r0*r1+(-2+3n+3n^2)*r1^2)/(-6-4n+24n^2+16n^3))
	fac = 1/sqrt( l*(l+1)*(-((1+l+l^2+2n*(1+n+n^2))*r0^2)+2(1+l+l^2+n-n^2)*r0*r1-(1+l+l^2+2n*(1+n+n^2))*r1^2)/((1+2n)*(r0-r1)) )
	# return (1-x)*(x+1)*jacobi(n-2,0,l+1/2, x) # - jacobi(n-1,0.0,0.0, x)
	return fac*(r1-r)*(r0-r)*jacobi(n-2,0,0, x) # - jacobi(n-1,0.0,0.0, x)
end


# @inline _nrange_p(b::Basis{Inviscid, SphericalShell},l) = 2:(b.N-l+1)
@inline _nrange_p(b::Basis{Inviscid, SphericalShell},l) = 2:(b.N-l+1)
@inline _nrange_t(b::Basis{Inviscid, SphericalShell},l) = 0:(b.N-l)

# @inline nrange_p_bc(b::Basis{Inviscid, SphericalShell},l) = 0:((b.N-l+1)÷2-1)
# @inline nrange_t_bc(b::Basis{Inviscid, SphericalShell},l) = nrange_t(b, l)

@inline function bcs_p(b::Basis{Inviscid, SphericalShell})
    # fs = (@inline((l, n) -> s(Basis{Inviscid, SphericalShell}, b.V, l, 0, n, b.V.r0)), 
	# 	  @inline((l, n) -> s(Basis{Inviscid, SphericalShell}, b.V, l, 0, n, b.V.r1)))
	fs = ()
    return fs
end

@inline function bcs_t(b::Basis{Inviscid, SphericalShell})
    fs = ()
	return fs
end


lpmax(b::Basis{Inviscid, SphericalShell}) = b.N
ltmax(b::Basis{Inviscid, SphericalShell}) = b.N


end