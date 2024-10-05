---
title: 'Limace.jl: A Julia package to compute hydromagnetic modes in spherical domains'
tags:
  - Julia
  - geophysics
  - fluid dynamics
  - geomagnetism
  - hydromagnetic modes
authors:
  - name: Felix Gerick
    orcid: 0000-0001-9924-0562
    affiliation: "1, 2"
affiliations:
 - name: National Centre for Space Studies, France
   index: 1
   ror: 04h1h0y33
 - name: Royal Observatory of Belgium, Belgium
   index: 2
   ror: 00hjks330
date: 2 October 2024
bibliography: paper.bib
---

# Summary

Hydromagnetic modes in spherical domains are relevant to the liquid cores of planets, moons or stars, as well as rotating fluid dynamics experiments.
These modes are solutions to the linearized rotating magnetohydrodynamic equations that govern electrically conducting fluids under rapid rotation. 
`Limace.jl` is a package written in the Julia programming language [@bezansonjulia2017], based on Galerkin projections of the governing equations onto 
trial vectors of the velocity and magnetic field. It aims to facilitate the calculation of modes in flexible setups with a high-level interface,
whilst remaining computationally performant enough to tackle relevant physical parameters.

# Statement of need

- Recent research interest in modelling these modes [@gerickfast2021, @trianacore2022, @luowaves2022a, @luowaves2022, @gerickinterannual2024b]
- No open-source framework for hydromagnetic modes in arbitrary background field geometry 
- To my knowledge, the only open source code to compute hydromagnetic modes in planetary cores is Kore [@trianaviscous2021].
- Several in-house closed source codes exist in the community
- Reimplementation requires substantial effort due to complexity of the spectral equations
- Generic and does not rely on symmetry assumptions.

# Theoretical background and implementation details

The modeled equations are based on the work of [@iversscalar2008] and [@gerickinterannual2024]. 

Integrals over the spherical surfaces are computed through the Adam-Gaunt and Elsasser variables [@jamesadams1973], which are calculated from Wigner symbols (available in Julia through [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl), based on [@johanssonfast2016]).
The remaining integration in radial direction is done using Gauss-Legendre quadratures available through [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).

To compute sparse solutions, a shift-invert spectral transform method is provided, based on the sparse LU factorization from `UMFPACK` [@davisalgorithm2004] and the partial Schur decomposition implemented in [ArnoldiMethod.jl](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl) [@StoppelsArnoldiMethod].

For postprocessing, `Limace.jl` uses a fast spherical harmonic transform implemented in the `SHTns` library [@schaefferefficient2013], and available in Julia through [SHTns.jl](https://github.com/fgerick/SHTns.jl).
It is used to transform the spectral coefficients to vector fields evaluated on a spatial grid.



# Basic Example - Malkus background magnetic field

In this example, we calculate the spectrum of modes when the background field is ``\mathbf{B}_0 = s\mathbf{e}_\phi``, following [@malkushydromagnetic1967]. 

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

λ = eigvals(Matrix(RHS))
```

In this simple configuration analytical solutions can be derived, relating the frequency of the hydromagnetic problem to the inertial mode frequencies.
We verify that our calculated mode spectrum contains some of the analytical solutions:

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

More examples are given in the documentation of `Limace.jl`.

# Acknowledgements

I have received funding from the European Research Council (ERC) GRACEFUL Synergy Grant No. 855677. 
This project has been funded by ESA in the framework of EO Science for Society, through contract 4000127193/19/NL/IA (SWARM + 4D Deep Earth: Core). 
I thank Phil Livermore for the key contributions in the theoretical development of the model.

# References
