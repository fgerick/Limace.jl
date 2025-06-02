# Bases

Several bases are readily implemented for full sphere geometries and satisfying certain boundary and orthogonality conditions. 
The most commonly encountered ones are outlined here and a small introduction into how to write a new basis is given at the end.

Each basis essentially consists of two scalar functions `s` and `t`, corresponding to the radial poloidal and toroidal scalars, respectively. 
They are introduced in the [Theoretical background](@ref) in detail as $S_{ln}$ and $T_{ln}$. 
The function arguments of `s` and `t` include the azimuthal degree `m` to keep it generic for possible future extensions.

Currently, all bases rely on a radial representation using Jacobi polynomials $J_n^{\alpha,\beta}(x)$ (implemented as [Limace.Poly.jacobi](@ref)).

All implemented bases are of type `Basis{T,Vol<:Volume}`, where the type `T` defines the type of the basis (e.g. `Inviscid`).
```@docs
Limace.Bases.Basis
```

## Inviscid velocity basis

For an inviscid fluid, we only require ``\mathbf{u}\cdot\mathbf{n}=0`` at the boundary. This is known as the non-penetration condition. For poloidal and toroidal decomposition, this boils down to the requirement that the poloidal scalar vanishes at the boundary. The toroidal scalar remains unconstrained and we have free slip.

The basis is orthonormal on the unit sphere $\mathcal{V}$, so that
```math
\int_\mathcal{V} \mathbf{u}_i^*\cdot\mathbf{u}_j\,\mathrm{d}V = \delta_{ij}.
```

```@autodocs
Modules = [Limace.InviscidBasis]
```


To define the inviscid basis, using a polynomial truncation `N = 10`, simply run
```@example bases
using Limace #hide
N = 10
u = Inviscid(N)
```

The basis `u` will include all azimuthal wave numbers `m`. To restrict to a single wave number (e.g. `m=1`), one can do

```@example bases
u = Inviscid(N; m=1)
```

## Viscous velocity basis 

For a viscous fluid, in addition to the non-penetration condition, we require no-slip, so that in total ``\mathbf{u} = \mathbf{0}`` at the boundary. [chen_optimal_2018](@citet) have written a basis that satisfies these conditions, and requires orthogonality w.r.t the vector Laplacian, so that
```math
\int_\mathcal{V} \mathbf{u}_i^*\cdot\boldsymbol{\nabla}^2\mathbf{u}_j\,\mathrm{d}V = \delta_{ij}.
```
This is desirable, as the resulting projections of all considered forces are banded (if only combined with another basis that satisfies the appropriate orthogonality).

```@autodocs
Modules = [Limace.ViscousBasis]
```



To define the viscous basis, using a polynomial truncation `N = 10`, simply run
```@example bases
u = Viscous(N)
```

The basis `u` will include all azimuthal wave numbers `m`. To restrict to a single wave number (e.g. `m=1`), one can do

```@example bases
u = Viscous(N; m=1)
```

## Insulating magnetic field basis

This basis can be used when the exterior is assumed to be a perfect insulator.
The magnetic field is then required to be continuous through the boundary, i.e.
```math
\left[\mathbf{B}\right]_{\delta\mathcal{V}} = 0,
```
where ``\left[x\right]_{\delta\mathcal{V}}`` denotes a jump across the boundary. 
In addition, the magnetic field is required to vanish at an infinite distance to the origin and to match a potential field in the exterior (so that ``\boldsymbol{\nabla}\times\mathbf{B} = \mathbf{0}`` in the exterior).
One can derive the required boundary condition at the outer boundary on the poloidal and toroidal scalars to be
```math
\begin{align*}
\left(\frac{\partial s_{lmn}}{\partial r} + \frac{l+1}{r}s_{lmn}\right)_{\partial\mathcal{V}} &= 0,\\
t_{lmn}\bigg|_{\partial\mathcal{V}} &= 0.
\end{align*}
```
for all $l,m,n$.

The basis is orthonormal w.r.t the vector Laplacian, so that
```math
\int_{\mathbb{R}^3} \mathbf{B}_i^*\cdot\boldsymbol{\nabla}^2\mathbf{B}_j\,\mathrm{d}V = \delta_{ij}.
```

```@autodocs
Modules = [Limace.InsulatingBasis]
```

To define the insulating basis, using a polynomial truncation `N = 10`, simply run
```@example bases
b = Insulating(N)
```

The basis `b` will include all azimuthal wave numbers `m`. To restrict to a single wave number (e.g. `m=1`), one can do

```@example bases
b = Insulating(N; m=1)
```

## Perfectly conducting magnetic field basis

This basis is appropriate for the magnetic field, when the exterior is a perfect conductor. In that case, the boundary condition is identical to the non-penetration condition for an inviscid flow. 

```math
\mathbf{B}\cdot\mathbf{n}\bigg|_{\delta \mathcal{V}} = 0,
```
which is equivalent to the condition

```math
s\bigg|_{\delta \mathcal{V}} = 0.
```

The basis is orthonormal on the unit sphere $\mathcal{V}$, so that
```math
\int_\mathcal{V} \mathbf{B}_i^*\cdot\mathbf{B}_j\,\mathrm{d}V = \delta_{ij}.
```

```@autodocs
Modules = [Limace.PerfectlyConductingBasis]
```


To define the perfectly conducting basis, using a polynomial truncation `N = 10`, simply run
```@example bases
b = PerfectlyConducting(N)
```

The basis `b` will include all azimuthal wave numbers `m`. To restrict to a single wave number (e.g. `m=1`), one can do

```@example bases
b = PerfectlyConducting(N; m=1)
```

## Unconstrained

This basis does not impose any boundary condition. It is a useful starting point for imposing custom boundary conditions explicitly in the linear operator.

The basis is orthonormal on the unit sphere $\mathcal{V}$, so that
```math
\int_\mathcal{V} \mathbf{u}_i^*\cdot\mathbf{u}_j\,\mathrm{d}V = \delta_{ij}.
```

```@autodocs
Modules = [Limace.UnconstrainedBasis]
```


To define the unconstrained basis, using a polynomial truncation `N = 10`, simply run
```@example bases
u = Unconstrained(N)
```

The basis `u` will include all azimuthal wave numbers `m`. To restrict to a single wave number (e.g. `m=1`), one can do

```@example bases
u = Unconstrained(N; m=1)
```


## Implementing a new basis

All bases are defined in a submodule, as parametric types `T` of a [Basis](@ref) struct and associated with a [BoundaryCondition](@ref).


We can follow the example of `Limace.InsulatingBasisNoBC`, which implements an insulating magnetic field basis, with the boundary condition explicitly imposed afterwards.

Define a new module in a `.jl` file in the `bases` folder, and import several components of `Limace.jl` that are needed:
```julia
module InsulatingBasisNoBC

using SparseArrays


using SparseArrays
using LinearAlgebra
using DocStringExtensions

using ..Bases
using ..UnconstrainedBasis, ..InviscidBasis, ..InsulatingBasis
using ..Utils
using ..Poly

using ..Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax, Sphere
import ..Bases: lpmax, ltmax, lmn_t, lmn_p, _nrange_p, _nrange_t, np, nt, t, s, nrange_p_bc, nrange_t_bc, bcs_p, bcs_t
using ..Poly: ∂
```

We define our new basis, with the right boundary condition. Here, we choose `BC=NoBC()`, as we impose the boundary condition explicitly afterwards.

```julia
export InsulatingNoBC

struct InsulatingNoBC end

InsulatingNoBC(N; kwargs...) = Basis{InsulatingNoBC,Sphere}(; N, BC=NoBC(), V=Sphere(), kwargs...)
```

We need to define a method for `s` and `t` for our basis. In this example, we just use the `Unconstrained` basis and add boundary conditions explicitly (i.e. they are not included in the basis elements).

```julia
s(::Type{Basis{InsulatingNoBC,Sphere}}, V::Volume, l,m,n,r)  = s(Basis{Unconstrained}, V, l,m,n,r) 
t(::Type{Basis{InsulatingNoBC,Sphere}}, V::Volume, l,m,n,r)  = t(Basis{Unconstrained}, V, l,m,n,r) 
```

A few additional functions need to be defined to determine the radial and spherical harmonic degrees.

```julia
@inline _nrange_p(b::Basis{InsulatingNoBC,Sphere}, l) = 0:((b.N-l+1)÷2)
@inline _nrange_t(b::Basis{InsulatingNoBC,Sphere}, l) = 0:((b.N-l)÷2)

@inline lpmax(b::Basis{InsulatingNoBC,Sphere}) = b.N
@inline ltmax(b::Basis{InsulatingNoBC,Sphere}) = b.N
```

If we have no boundary condition imposed in the poloidal and toroidal scalars (`NoBC()`), we need to define the explicit boundary condition functions for the toroidal (`bcs_t`) and poloidal (`bcs_p`) scalars.

```julia
@inline function bcs_t(b::Basis{InsulatingNoBC,Sphere})
    fs = (@inline((l, n) -> t(Basis{InsulatingNoBC,Sphere}, b.V, l, 0, n, b.V.r1)),)
    return fs
end

@inline function bcs_p(b::Basis{InsulatingNoBC,Sphere})
    fs = (@inline((l, n) -> ∂(r -> s(Basis{InsulatingNoBC,Sphere}, b.V, l, 0, n, r), b.V.r1) + (l + 1) * s(Basis{InsulatingNoBC,Sphere}, b.V, l, 0, n, b.V.r1)),)
    return fs
end

end #module
```

Here, these are given as evaluations at the surface, `b.V.r1 = 1.0`. We use automatic differentiation to compute the derivatives in `r` using `ForwardDiff.jl`.

## Boundary conditions

The boundary condition is indicated in the basis, e.g.

```@example bases
u = Inviscid(N)
u.BC
```

Currently, only [Limace.Bases.NoBC](@ref) has a significant difference to the other boundary conditions.  

```@docs
Limace.Bases.NoBC
```

When a basis has a `NoBC` boundary condition, it considers that the basis elements do not contain a boundary condition themselves.
During the assembly of the linear operators (as introduced in [Implementing a new basis](@ref)), explicit lines for the boundary conditions will be left empty.
The boundary conditions then need to be added separately, as outlined in [Explicit boundary condition operator](@ref).

All boundary conditions are a subtype of
```@docs
Limace.Bases.BoundaryCondition
```

The other "placeholder" boundary conditions, currently implemented, are
```@docs
Limace.Bases.InsulatingBC
Limace.Bases.InviscidBC
Limace.Bases.NoSlipBC
Limace.Bases.PerfectlyConductingBC
```

## Basis elements

Single elements of a basis are wrapped into a `BasisElement` structure
```@docs
Limace.Bases.BasisElement
```

