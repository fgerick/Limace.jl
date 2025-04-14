```@meta
CurrentModule = Limace
```

# Limace.jl

Documentation for [Limace.jl](https://github.com/fgerick/Limace.jl), a linear model based on Galerkin projections to solve for hydromagnetic modes in rotating spheres.

## Installation
To get started, you will need a working Julia environment, preferrably `>=v1.10`. Then, from the REPL, you can run 

```julia
import Pkg; Pkg.add("https://github.com/fgerick/Limace.jl.git")
```

This should install most of what you need to compute the basics. 
To visualize the solutions you will need to install a plotting package of your choice.
The examples here use [Makie.jl](https://docs.makie.org/dev/) and [GeoMakie.jl](https://github.com/MakieOrg/GeoMakie.jl). 

## Getting started

To solve for modes, four steps are needed:

1) Definition/Choice of appropriate bases.
2) Chose forcings, background state and parameters.
3) Assembly of Galerkin projection matrices.
4) Compute solution(s) of (generalized) eigen problem.

These four steps, are accessible through a high-level interface, as presented in the [Quickstart](@ref).
The details are introduced in a more comprehensive way through the Examples, and in the relevant sections.

Some [Theoretical background](@ref) is given, as well as some more detailed API information on the [Bases](@ref) and the implemented [Forces](@ref).