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

To solve for modes, three steps are needed:

1) Definition/Choice of appropriate bases.
2) Assembly of Galerkin projection matrices.
3) Compute solution(s) of (generalized) eigen problem.

These three steps, with some post-processing, are introduced best through the Examples.

**But I just want to copy-paste code to see if it works**

Calculate hydromagnetic modes for the Malkus field:

```julia
using Limace, LinearAlgebra

N = 6
u = Inviscid(N)
b = PerfectlyConducting(N) # == Inviscid(N)

B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15)) # corresponds to B_0 = s e_phi	
Le = 1e-2

RHSc = Limace.coriolis(u)
RHSl = Limace.lorentz(u, b, B0)
RHSi = Limace.induction(b,u,B0)
RHSd = spzeros(length(b),length(b)) #empty (no Ohmic diffusion)
RHS = [RHSc/Le RHSl
	   RHSi RHSd];

λ, x = eigen(Matrix(RHS))
```

Some [Theoretical background](@ref) is given, as well as some more detailed API information on the [Bases](@ref) and the implemented [Forces](@ref).