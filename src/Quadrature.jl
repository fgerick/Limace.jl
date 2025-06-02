module Quadrature

using FastGaussQuadrature
using ..Bases: Volume

export rquad, ∫dr

# ∫₀¹ f(r) r² dr
function rquad(nr) 
	# _r, _wr = gaussjacobi(nr,0.0,2.0) 
	_r, _wr = gausslegendre(nr)
	_r = (_r .+ 1)/2
	_wr /= 2
	# _wr /= 8
	return _r, _wr
end

# ∫ₐᵇ f(r) r² dr
function rquad(nr,a,b) 
	# _r, _wr = gaussjacobi(nr,0.0,2.0) 
	_r, _wr = gausslegendre(nr)
	_r = ((b-a)*_r .+ b .+ a)/2
	_wr *= (b-a)/2
	# _wr /= 8
	return _r, _wr
end

rquad(nr, V::Volume) = rquad(nr, V.r0, V.r1)

"""
    rquad(nr, a, b)
    rquad(nr, V::Volume)

Computes `nr` Gauss-Legendre quadrature points and weights on a radial grid. 

```julia
r0 = 0.3
r1 = 1.0

r, wr = Limace.Quadrature.rquad(50, r0, r1)

V = SphericalShell(r0,r1)

r, wr = Limace.Quadrature.rquad(50, V)
```

"""
function rquad end

"""
    ∫dr(f, r::Vector{T}, wr::Vector{T})::Complex{T} where T

Computes the integral of a function `f`
```math
\\int_{r_0}^{r_1} r^2 f(r)\\,\\mathrm{d}r,
```

with ``r_0<r<r_1``, using quadrature points `r` and weights `wr` computed by [Limace.Quadrature.rquad](@ref).
This function returns a complex value by default for easier type stability.

Type `\\int<TAB>` to write ∫.
"""
function ∫dr(f::F,r::Vector{T},wr::Vector{T})::Complex{T} where {F,T}
	out = zero(ComplexF64)
	for (r,w) in zip(r,wr)
		out+=f(r)*r^2*w
	end
	return out
end



end #module
