# Quickstart through high-level interface

`Limace.jl` provides a high-level interface to setup and solve the linear problem. 
The four main steps to solve for hydromagnetic modes are outlined here:

## Definition/Choice of appropriate bases.

Chose a polynomial truncation, e.g. `N=10`, and an appropriate bases (see [Bases](@ref) for implemented bases) for the velocity (and magnetic field).
Let us use the [`Inviscid`](@ref "Inviscid velocity basis") basis for the velocity and the [`PerfectlyConducting`](@ref "Perfectly conducting magnetic field basis") for the magnetic field.

```julia
N = 10
u = Inviscid(N)
b = PerfectlyConducting(N)
bases = [u,b]
```

## Chose forcings, background state and parameters.

We can chose a background magnetic field (or flow) and an appropriate non-dimensional parameter (e.g. Lehnert number)
```julia
B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15)) # corresponds to B_0 = s e_z
Le = 1e-2 #Lehnert number
```

The forcings considered in the problem are collected in a list to be included
```julia
forcings = [Limace.Inertial(u), 
	Limace.Inertial(b), 
	Limace.Coriolis(u, 1/Le), 
	Limace.Lorentz(u,b,B0), 
	Limace.InductionB0(b,u,B0)]
```

We can combine the `bases` and `forcings` in a [`Limace.LimaceProblem`](@ref).
```julia
problem = LimaceProblem(bases, forcings)
```

```@docs
Limace.LimaceProblem
```
## Assembly of Galerkin projection matrices.

The assembly of the linear operators within the `problem` is then simply:
```julia
Limace.assemble!(problem)
```
We can use a keyword `threads=true` to assemble the matrices using multithreading on a single computer (see [`Limace.assemble!`](@ref)).
```@docs
Limace.assemble!
```

## Compute solution(s) of (generalized) eigen problem.

Solving the eigen problem can be done through a high-level function
```julia
Limace.solve!(problem)
```

We can then access the eigenvalues `λ`, and eigenvectors `x` within `problem.sol`:
```julia
λ, x = problem.sol.values, problem.sol.vectors
```

```@docs
Limace.solve!
```

More details on the choices of method and the underlying methodology is outlined in the section on [Solving the eigenvalue problem](@ref "Solving the eigenvalue problem").


From here, we can plot the solutions or process the spectrum of modes. See the examples and [Postprocessing](@ref) for more details.
