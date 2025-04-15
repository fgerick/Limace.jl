# Solving the eigenvalue problem

The generalized eigenvalue problem is given as
```math
\lambda \mathbf{B}\mathbf{x} = \mathbf{A}\mathbf{x},
```
where ``\lambda`` is the eigenvalue and ``\mathbf{x}`` the corresponding (right) eigenvector.
The linear operators (or matrices) ``\mathbf{A}`` and ``\mathbf{B}`` contain the Galerkin-projected operators of the system of equations.
In the context of `Limace.jl`, these matrices are often referred to as `LHS` (_"left-hand-side"_) for ``\mathbf{B}``, and `RHS` (_"right-hand-side"_) for ``\mathbf{A}``.
When ``\mathbf{B}`` is the unit matrix, the eigenvalue problem reduces to the so-called standard eigenvalue problem. This is for example the case when ``\mathbf{B}`` only contains the inner products of orthonormal bases.

Solving the generalized eigenvalue problem is not a trivial task for matrices that are large (already at a dimension ``\sim 10`` this becomes difficult analytically).
For different sizes of ``\mathbf{A}`` and ``\mathbf{B}``, we can/have to rely on different numerical methods to solve for eigensolutions ``(\lambda,\mathbf{x})``.

## Dense vs. sparse solutions

When the matrices ``\mathbf{A}`` and ``\mathbf{B}`` are small enough, it is possible to represent them as dense matrices and to compute the full (_dense_) spectrum of eigenvalues and eigenvectors. 
_"Small enough"_ depends on the computer you are using, but for matrix sizes larger than ``10^4``, solving for all eigenvalues (and eigenvectors) is no longer feasible even on the average workstation. 
Such matrix sizes are quickly approached for truncation degrees ``N\geq 40`` in `Limace.jl`, when no symmetry assumptions can be used (e.g. we cannot consider azimuthal degrees independently due to non-axisymmetry).

Luckily, there are iterative methods to solve only for a subset of eigenvalues of the large sparse matrix system. 
Several different methods exist, with one widely used method being the Arnoldi method, which relies on the orthonormalization of Krylov subspaces. 
Through these methods we can compute a few eigenvalues at the extremes of the spectrum (e.g. the largest eigenvalue in absolute value), and also eigenvalues close to a target by shifting the spectrum.

How to access these different eigen solvers is outlined in the following.

## Dense eigenvalue solver

The dense eigensolver routines provided through [LAPACK](https://www.netlib.org/lapack/) are available in Julia through the standard library `LinearAlgebra`:
```julia
using LinearAlgebra
```

Let us consider the torsional mode problem (see examples), which assembles the two matrices `LHS` and `RHS`, so that
```julia
using Limace

N = 40
Le = 5e-4
Lu = 1/Le
B0 = [BasisElement(Basis{Insulating}, Poloidal, (1,0,1),0.3), BasisElement(Basis{Insulating}, Poloidal, (2,0,1),0.7)]

u = Inviscid(N; m=0)
b = Insulating(N; m=0)
bases = [u,b]
forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u,1/Le), Limace.Lorentz(u,b,B0), Limace.InductionB0(b,u,B0), Limace.Diffusion(b,1/Lu)]
problem = LimaceProblem(bases,forcings)
Limace.assemble!(problem; threads=true)
```

We can then solve the dense [`LimaceProblem`](@ref) by using
```julia
Limace.solve!(problem)
```
which uses the default keyword argument `method=:dense`.

The results are stored in `problem.sol`, so that
```julia
λ, x = problem.sol.values, problem.sol.vectors
```

The function [`Limace.solve!`](@ref) essentially converts the sparse matrices `LHS` and `RHS` of size `2500x2500` to dense matrices, using
```julia
(; LHS, RHS) = problem # syntactic sugar for LHS,RHS = problem.LHS, problem.RHS
LHSd = Matrix(LHS)
RHSd = Matrix(RHS)
```

and solves for the eigensolutions
```julia
λ, x = eigen(RHSd, LHSd)
```

to get all eigenvalues `λ` and all eigenvectors `x`, so that

```julia
i=1
λ[i]*LHS*x[:,i] ≈ RHS*x[:,i]
```
for all `i`.

If one requires only the eigenvalues, we can save computational effort by using
```julia
λ = eigvals(RHSd, LHSd)
```


If `LHS` is a unit matrix, we can solve the problem as
```julia
λ, x = eigen(RHSd)
```

## Sparse eigenvalue solver

A native Julia implementation of the Arnoldi method (with Krylov-Schur restarts) is implemented in [ArnoldiMethod.jl](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl).
In the `EigenSolve` submodule of `Limace.jl`, an implementation of the shift-invert method is available through the [`Limace.EigenSolve.eigstarget`](@ref) function.
The shift-invert method shifts the spectrum around a given target eigenvalue ``\sigma`` and inverts the operator on the left-hand-side, so that
```math
\frac{1}{\lambda-\sigma}\mathbf{x} = (\mathbf{A}-\sigma\mathbf{B})^{-1}\mathbf{B}\mathbf{x}.
```
Which is a standard eigenvalue problem
```math
\hat{\lambda}\mathbf{x} = \mathbf{C}\mathbf{x},
```
with ``\hat{\lambda} = (\lambda-\sigma)^{-1}``, and ``\mathbf{C} = (\mathbf{A}-\sigma\mathbf{B})^{-1}\mathbf{B}``. 
When computing the largest-amplitude eigenvalues to this problem, we find the eigenvalues closest in amplitude to the target eigenvalue ``\sigma``.


```@docs
Limace.EigenSolve.eigstarget
```

Using the same matrices assembled previously, we can calculate solutions close to a target ``\sigma = \mathrm{i}`` (frequency ``\omega=1``).
```julia
target = 1.0im

λ, x = Limace.EigenSolve.eigstarget(RHS, LHS, target; nev=5);
```

this is equivalent to using the high-level interface

```julia
Limace.solve!(problem; method=:sparse, target=1.0im, nev=5)
```

The standard inverse method is also provided through [`Limace.EigenSolve.eigs`](@ref), in order to compute only the extremal eigenvalues (no shift required).
This is approach is not possible when ``\mathbf{B}`` is singular (e.g. when explicit boundary conditions are imposed, ``\mathbf{B}`` contains zero lines and is therefore singular).

```@docs
Limace.EigenSolve.eigs
```

For example we can compute the largest eigensolution using the high-level interface
```julia
Limace.solve!(problem; method=:sparse, nev=1)
```

or 

```julia
λ, x = Limace.EigenSolve.eigs(RHS, LHS; nev=1);
```

These sparse approaches are computationally much more efficient and feasible than the dense approach for much larger truncation degrees `N` (as we also only compute a few eigensolutions and not the complete spectrum).
Still, computational effort increases as the number of requested eigensolutions `nev` is raised, or the matrix size is increased. 


