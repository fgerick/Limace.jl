# Forces

Different forces included in `Limace.jl`. 

So far, there is no heat equation in `Limace.jl`, so the forces include only the mass-term (or `Limace.inertial`), the Coriolis force (`Limace.coriolis`), the Lorentz force (`Limace.lorentz`) and diffusion (`Limace.diffusion`). 
In the induction equation we have the magnetic induction term $\nabla\times\mathbf{u}\times\mathbf{B}$ (`Limace.induction`). Some details for each of these forcings is given on the respective pages.

## Inertial

```math
\int \mathbf{u}_i^* \cdot \mathbf{u}_j\,\mathrm{d}V,
```

```@autodocs
Modules = [Limace]
Pages = ["forces/inertial.jl"]
```

## Coriolis

The necessary functions to compute Galerkin projections on the Coriolis term

```math
\int \mathbf{u}_i^* \cdot 2\mathbf{e}_z\times\mathbf{u}_j\,\mathrm{d}V,
```

where ``\mathbf{u}_{i,j}`` are either poloidal or toroidal basis vectors.

```@autodocs
Modules = [Limace]
Pages = ["forces/coriolis.jl"]
```

## Diffusion 

```math
\int \mathbf{u}_i^* \cdot \boldsymbol{\nabla}^2\mathbf{u}_j\,\mathrm{d}V
```

```@autodocs
Modules = [Limace]
Pages = ["forces/diffusion.jl"]
```

## Induction

```math
\int \mathbf{b}_i^* \cdot \boldsymbol{\nabla}\times\left(\mathbf{u}_j\times\mathbf{B}_k\right)\,\mathrm{d}V
```

```@autodocs
Modules = [Limace]
Pages = ["forces/induction.jl"]
```

## Lorentz

```math
\int \mathbf{b}_i^* \cdot \left(\boldsymbol{\nabla}\times\mathbf{b}_j\times\mathbf{b}_k\right)\,\mathrm{d}V
```

```@autodocs
Modules = [Limace]
Pages = ["forces/lorentz.jl"]
```


## Explicit boundary condition operator

```@autodocs
Modules = [Limace]
Pages = ["forces/bc.jl"]
```