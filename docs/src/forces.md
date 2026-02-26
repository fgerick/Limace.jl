# Forces

Different forces (or forcings, linear operators in general) are included in `Limace.jl`.
For the high-level interface, they are implemented as abstract types `Limace.Forcing{1}` or `Limace.Forcing{2}`, for on one or two bases dependencies respectively.

Currently, the implemented forces are:
- [Inertial](#Inertial)
- [Advection](#Advection)
- [Coriolis](#Coriolis)
- [Diffusion](#Diffusion)
- [Induction](#Induction)
- [Lorentz](#Lorentz)

Some details for each of these forces is given in the following.

## Inertial

```math
\int \mathbf{u}_i^* \cdot \mathbf{u}_j\,\mathrm{d}V,
```

The linear operator is simply the inner product of the basis elements. For the momentum equation, this is essentially the projection of the inertial term ``\partial_t \mathbf{u} \rightarrow \lambda \mathbf{u}``, hence the name.

```@docs 
Limace.Inertial
Limace.inertial
```

## Advection

```math
\int \mathbf{u}_i^* \cdot \left(\left(\boldsymbol{\nabla}\times\mathbf{u}_j\right)\times\mathbf{U}_0 + \left(\boldsymbol{\nabla}\times\mathbf{U}_0\right)\times\mathbf{u}_j \right)\,\mathrm{d}V
```

This equivalent to the projection of the operator $\left(\mathbf{u}\cdot\boldsymbol{\nabla}\right)\mathbf{U}_0+\left(\mathbf{U}_0\cdot\boldsymbol{\nabla}\right)\mathbf{u}$, when $\boldsymbol{\nabla}\cdot\mathbf{u} = \boldsymbol{\nabla}\cdot\mathbf{U}_0 = 0$.
```@docs
Limace.Advection
Limace.advection
```

## Buoyancy

```math
\int \mathbf{u}_i^* \cdot \left(T_j\mathbf{r} \right)\,\mathrm{d}V
```

```@docs
Limace.Buoyancy
Limace.buoyancy
```

## Coriolis


```math
\int \mathbf{u}_i^* \cdot 2\mathbf{e}_z\times\mathbf{u}_j\,\mathrm{d}V,
```

```@docs
Limace.Coriolis
Limace.coriolis
```

## Diffusion 

```math
\int \mathbf{u}_i^* \cdot \boldsymbol{\nabla}^2\mathbf{u}_j\,\mathrm{d}V
```

```@docs
Limace.Diffusion
Limace.diffusion
```

## Induction

```math
\int \mathbf{b}_i^* \cdot \boldsymbol{\nabla}\times\left(\mathbf{u}_j\times\mathbf{B}_0\right)\,\mathrm{d}V
```

```@docs
Limace.InductionB0
```

```math
\int \mathbf{b}_i^* \cdot \boldsymbol{\nabla}\times\left(\mathbf{U}_0\times\mathbf{b}_j\right)\,\mathrm{d}V
```

```@docs
Limace.InductionU0
```

```@docs
Limace.induction
```



## Lorentz

```math
\int \mathbf{u}_i^* \cdot \left(\left(\boldsymbol{\nabla}\times\mathbf{b}_j\right)\times\mathbf{B}_0+\left(\boldsymbol{\nabla}\times\mathbf{B}_0\right)\times\mathbf{b}_j\right)\,\mathrm{d}V
```

```@docs
Limace.Lorentz
Limace.lorentz
```

## Scalar advection

```math
\int T_i^* \left(\mathbf{u}_j\boldsymbol{\cdot}\nabla T_0\right)\,\mathrm{d}V
```

```@docs
Limace.ScalarAdvectionT0
```

```math
\int T_i^* \left(\mathbf{U}_0\boldsymbol{\cdot}\nabla T_j\right)\,\mathrm{d}V
```

```@docs
Limace.ScalarAdvectionU0
```
```@docs
Limace.scalaradvection
```

## Explicit boundary condition operator

```@autodocs
Modules = [Limace]
Pages = ["forces/bc.jl"]
```

## Assembly

All forcings are constructed without being assembled (i.e. the projections of each forcing onto the relevant basis are not done at the definition of the forcing).

The projections are done explicitly through `Limace.assemble!`, that saves the assembled sparse matrix in the forcing structure (`f`).

```@example forces
using Limace #hide

N = 10
u = Inviscid(N)
b = Insulating(N)
B0 = BasisElement(b, Poloidal, (1,0,1))
f = Limace.Lorentz(u,b,B0)
Limace.assemble!(f)
```

This can also be done in parallel, by using the keyword argument `threads=true`.

```@example forces
Limace.assemble!(f; threads=true)
nothing #hide
```

When the forcings are wrapped into a [LimaceProblem](@ref), it is not necessary to call `Limace.assembly!` on each forcing by hand, but one can simply call `Limace.assembly!(problem)`, as outlined in the [Quickstart](@ref "Quickstart through high-level interface").
