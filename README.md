# Limace.jl
<!-- _**Li**near **m**odes **a**t the **c**ore of **E**arth_ -->
_**L**inear **I**nertial **MA**gneto **C**oriolis **E**igenmodes_
<!-- _**L**inear **I**nertial **M**agneto **A**rchimedes **C**oriolis **E**igenmodes_ -->

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://fgerick.github.io/Limace.jl/stable/) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://fgerick.github.io/Limace.jl/dev/) [![Build Status](https://github.com/fgerick/Limace.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/fgerick/Limace.jl/actions/workflows/CI.yml?query=branch%3Amain) [![Coverage](https://codecov.io/gh/fgerick/Limace.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/fgerick/Limace.jl)

<img src="limace_logo.jpg" width="40%">

Solve for linear eigenmodes of the rotating magnetohydrodynamics equations in a sphere.

## Installation
Simply run
```julia
import Pkg; Pkg.add("https://github.com/fgerick/Limace.jl.git")
```

## Examples 

### Inviscid inertial modes

For the inviscid inertial modes in a full sphere:
```julia
using Limace, LinearAlgebra

N = 10 #polynomial truncation
basis = Inviscid(N) #create Galerkin basis
A = Limace.coriolis(basis) #assemble Coriolis operator matrix (sparse)
λ = eigvals(Matrix(A)) #solve for eigenvalues
```
Note, the `Inviscid` basis is orthonormal, so that we do not need to calculate an operator associated with inertia.

We can check against analytical values available from the literature ([Zhang et al. 2001](https://doi.org/10.1017/S0022112001004049)):
```julia
function zhang(m, N) 
	sm = sign(m)
	m = abs(m)
	return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
end

any(λ .≈ zhang(1, 1)) #true
any(λ .≈ zhang(2, 1)) #true
any(λ .≈ zhang(3, 1)) #true

```

### Malkus modes

We create our two Bases for the flow and the magnetic field.
```julia

N = 6
u = Inviscid(N)
b = PerfectlyConducting(N) # == Inviscid(N)
```
Without specifying the azimuthal wave number $m$, e.g. 
```julia
u = Inviscid(N; m=1)
```
all $m \in [-l,l]$ with $l \in [1,N]$ are included. In the case of the Malkus field, this is not necessary, but in general (when $\mathbf{B}_0$ consist not only of $m=0$ components) we couple all $m$.

The background magnetic field $\mathbf{B}_0 = s \mathbf{e}_z$ is defined and we choose our characteristic time scale as the Alfvén time, so that our nondimensional parameter is the Lehnert number.
```julia
B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15)) # corresponds to B_0 = s e_z
Le = 1e-2
```
Then, we compute our projection operators for the Coriolis force, the Lorentz force and the induction term. In the ideal limit here, no diffusive term is included.
```julia
RHSc = Limace.coriolis(u)
RHSl = Limace.lorentz(u, b, B0)
RHSi = Limace.induction(b,u,B0)
RHSd = spzeros(length(b),length(u)) #empty (no Ohmic diffusion)
RHS = [RHSc/Le RHSl
	   RHSi RHSd];
```

Again, the projections on $\partial_t \mathbf{u}$ and $\partial_t \mathbf{b}$ are unit matrices, so that the eigenvalues are simply computed as
```julia
λ = eigvals(Matrix(RHS))
```
We can again compare to the analytical solutions:
```julia
function zhang(m, N) 
	sm = sign(m)
	m = abs(m)
	return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
end

# Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
slow(m, N, Le, λ = imag(zhang(m, N))) = im * λ / 2Le * (1 - √(1 + 4Le^2 * m * (m - λ) / λ^2))
fast(m, N, Le, λ = imag(zhang(m, N))) = im * λ / 2Le * (1 + √(1 + 4Le^2 * m * (m - λ) / λ^2))


for m = vcat(-(N-1):-1, 1:(N-1))
	@show any(isapprox(slow(m,1,Le)),λ)
	@show any(isapprox(fast(m, 1, Le)),λ)
end
```

More examples are in the [documentation](https://fgerick.github.io/Limace.jl/stable/) and the `test/modes.jl` file.

## Citation
