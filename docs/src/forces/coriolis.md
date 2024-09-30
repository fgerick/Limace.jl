# Coriolis force

The necessary functions to compute Galerkin projections on the Coriolis term

```math
\int \mathbf{u}_i \cdot 2\mathbf{e}_z\times\mathbf{u}_j\,\mathrm{d}V,
```

where ``\mathbf{u}_{i,j}`` are either poloidal or toroidal basis vectors.

```@autodocs
Modules = [Limace]
Pages = ["forces/coriolis.jl"]
```