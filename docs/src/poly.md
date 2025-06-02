
# Poly

Polynomial functions and conveniences that are needed for the spectral or spatial discretization.

## Polynomials

```@docs
Limace.Poly.jacobi
Limace.Poly.ylm
```

## Integral functions

Integrals over the spherical surfaces are reduced to two specific integrals, known as the Adam-Gaunt and Elsasser integrals [james_adams_1973](@citep).
```@docs
Limace.Poly.adamgaunt
Limace.Poly.elsasser
```

The radial part of the projection of poloidal and toroidal vectors can be described by two scalar functions, [Limace.Poly.inners](@ref) and [Limace.Poly.innert](@ref).
```@docs
Limace.Poly.inners
Limace.Poly.innert
```

These functions are used for the discrete integration in radius. 
The derivatives are computed using automatic differentiation (based on [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl)).

```@docs
Limace.Quadrature.rquad
```

The radial integrals of the projections are all integrated using
```@docs
Limace.Quadrature.∫dr
```

```@example poly
using Limace #hide
using Limace.Quadrature: rquad, ∫dr
r,wr = rquad(10, 0.0, 1.0)
f = r->r^2
∫dr(f,r,wr) ≈ 1/5	
```

## Miscellaneous functions

The functions listed here are mostly for internal convenience. 
Nevertheless, they may be useful for further implementations.

```@docs
Limace.Poly.p
Limace.Poly._∂ll
Limace.Poly.D
Limace.Poly.dylmdθ
Limace.Poly.dylmdϕ
```
