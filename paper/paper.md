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

The computation of hydromagnetic modes is under active scientific investigation, in particular their application in the study of the Earth's liquid core has recently gotten new attention [@gerickfast2021; @trianacore2022; @luowaves2022a; @luowaves2022; @gerickinterannual2024].
To model these waves, substantial work has to be put into the accurate implementation of the spectral equations that govern these modes.
This is true also for very basic and idealized examples, resulting in a substantial entry barrier for scientists to model these modes.
`Limace.jl` tries to lower this entry barrier, by providing an open source model with a very simple high-level API and modern online documentation with practical examples.

Despite having a high-level interface, `Limace.jl` can be used to solve complex and more geophysically relevant problem setups.
One of the only open-source models to compute hydromagnetic modes in planetary cores is Kore [@trianaviscous2021], 
a spectral code based on ultraspherical polynomials in axisymmetric setups written in Python.
A unique feature of `Limace.jl` is the support of complex background magnetic fields and flows over which the modes evolve.
The code has been developed from the beginning to leave assumptions of symmetry up to the user.
The model code base is tested against mode solutions from the scientific literature to ensure its accuracy.

# Theoretical background and implementation details

In non-linear form, the momentum equation of the incompressible fluid and the induction equation, are written as
$$
\begin{aligned}
	\frac{\partial \mathbf{U}}{\partial t} + \left(\boldsymbol{\nabla}\times\mathbf{U}\right)\times\mathbf{U} + 2\boldsymbol{\Omega}\times\mathbf{U} &= -\frac{1}{\rho}\nabla P + \frac{1}{\rho\mu_0} \left(\boldsymbol{\nabla}\times\mathbf{B}\right)\times\mathbf{B} + \nu\boldsymbol{\nabla}^2\mathbf{U} + \mathbf{F}, \\
	\frac{\partial \mathbf{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\mathbf{U}\times\mathbf{B}\right) + \eta \boldsymbol{\nabla}^2\mathbf{B},
\end{aligned}
$$
with $\mathbf{U}$ the velocity, $\mathbf{B}$ the magnetic field, $\boldsymbol{\Omega}$ the rotation axis, $\rho$ the fluid density, $P$ the reduced hydrodynamic pressure, $\mu_0$ the magnetic permeability of free space, $\nu$ the kinematic viscosity, $\mathbf{F}$ some additional body force, and $\eta$ the magnetic diffusivitiy.

In order to compute modal solutions, the velocity, magnetic and pressure fields a linearized, so that 
$$
\begin{aligned}
	\mathbf{U}(\mathbf{r},t) & = \mathbf{U}_0(\mathbf{r})+ \mathbf{u}(\mathbf{r}) e^{\lambda t}, \\
	\mathbf{B}(\mathbf{r},t) & = \mathbf{B}_0(\mathbf{r})+ \mathbf{b}(\mathbf{r}) e^{\lambda t}, \\
	P(\mathbf{r},t)   & = P_0(\mathbf{r})+ p(\mathbf{r}) e^{\lambda t}.
\end{aligned}
$$
with $\lambda=-\sigma+\mathrm{i}\omega$, with $\sigma$ the damping rate and $\omega$ the frequency of the oscillatory perturbation to the steady background.
Removing the steady part and neglecting higher order terms, the linearized momentum and induction equations then read
$$
\begin{aligned}
	\lambda\mathbf{u} =& -\left(\boldsymbol{\nabla}\times\mathbf{u}\right)\times\mathbf{U}_0- \left(\boldsymbol{\nabla}\times\mathbf{U}_0\right)\times\mathbf{u} -2\Omega\mathbf{e}_z\times\mathbf{u} - \frac{1}{\rho}\nabla p\\
	 &+ \frac{1}{\rho\mu_0}\left(\left(\boldsymbol{\nabla}\times\mathbf{b}\right)\times\mathbf{B}_0+\left(\boldsymbol{\nabla}\times\mathbf{B}_0\right)\times\mathbf{b}\right) + \nu \boldsymbol{\nabla}^2\mathbf{u},\nonumber\\
	\lambda\mathbf{b} =& \boldsymbol{\nabla}\times\left(\mathbf{U}_0\times\mathbf{b}\right) + \boldsymbol{\nabla}\times\left(\mathbf{u}\times\mathbf{B}_0\right) + \eta \boldsymbol{\nabla}^2\mathbf{b}.
\end{aligned}
$$

These equations are then projected onto trial vectors $\boldsymbol{\xi}_i$, so that
$$
f_{ij} = \int \boldsymbol{\xi}_i \cdot \mathbf{f}\left(\mathbf{u}_j,\mathbf{b}_j, \mathbf{U}_0, \mathbf{B}_0\right)\,\mathrm{d}V,
$$
where $\mathbf{f}$ is any of the terms in the momentum and induction equation and $\boldsymbol{\xi}_i = [\mathbf{u}_i, \mathbf{b}_i]$.
Due to the divergence free condition on the velocity and magnetic field, i.e. the flow is incompressible and no magnetic monopoles, 
it is convenient to decompose the fields into poloidal and toroidal components.
$$
\begin{aligned}
    \mathbf{u} &= \sum_i \alpha_i\mathbf{u}_i = \sum_{l,m,n} \alpha^P_{lmn}\mathbf{P}_{lmn} + \sum_{l,m,n} \alpha^Q_{lmn}\mathbf{Q}_{lmn},\\
    \mathbf{b} &= \sum_i \beta_i\mathbf{b}_i =  \sum_{l,m,n} \beta^S_{lmn}\mathbf{S}_{lmn} + \sum_{l,m,n} \beta^T_{lmn}\mathbf{T}_{lmn},
\end{aligned}
$$
with the respective poloidal and toroidal basis vectors
$$
\begin{aligned}
	\left[\mathbf{P},\mathbf{S}\right]_{lmn} & = \boldsymbol{\nabla}\times\boldsymbol{\nabla}\times \left[P,S\right]_{ln}(r)Y_l^m(\theta,\phi)\mathbf{r}, \\
	\left[\mathbf{Q},\mathbf{T}\right]_{lmn} & = \boldsymbol{\nabla}\times \left[Q,T\right]_{ln}(r) Y_l^m(\theta,\phi)\mathbf{r}.
\end{aligned}
$$
Here, $Y_l^m(\theta,\phi)$ is the (fully normalized) spherical harmonic of degree $l$ and order $m$.
The boundary conditions (or regularity condition at $r=0$) are imposed on the scalar functions $P,S,Q,T$.

We need to consider all combinations of poloidal and toroidal vector combinations in the projection of the forces.
This leads to several long coupling terms, especially for the Lorentz force and induction term. 
The integrals of these coupling terms over the spherical surfaces are computed through the Adam-Gaunt and Elsasser variables [@jamesadams1973], which are calculated from Wigner symbols (available in Julia through [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl), based on @johanssonfast2016).
The remaining integration in radial direction is done using Gauss-Legendre quadratures available through [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).
The exact modeled equations are outlined in @gerickinterannual2024, based on the work of @iversscalar2008.  

From the projected equations, the problem reduces to a generalized eigen problem
$$
\lambda \mathbf{A}\mathbf{x} = \mathbf{B}\mathbf{x},
$$
that is solved numerically. The matrix $\mathbf{B}$ is generally not symmetric/Hermitian, but $\mathbf{A}$ can be the unit matrix, symmetric tridiagonal or symmetric, depending on the chosen bases.

For small problem sizes, the eigen problem can be solved using dense methods, e.g. using the standard library function `eigen`.
To compute few eigen solutions of the sparse system, a shift-invert spectral transform method is provided, based on the sparse LU factorization from `UMFPACK` [@davisalgorithm2004] and the partial Schur decomposition implemented in [ArnoldiMethod.jl](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl) [@StoppelsArnoldiMethod].

For postprocessing, `Limace.jl` uses a fast spherical harmonic transform implemented in the [SHTns](https://bitbucket.org/nschaeff/shtns) library [@schaefferefficient2013], and available in Julia through [SHTns.jl](https://github.com/fgerick/SHTns.jl).
It is used to transform the spectral coefficients to vector fields evaluated on a spatial grid.



# Basic Example - Malkus background magnetic field

In this example, we calculate the spectrum of modes when the background field is $\mathbf{B}_0 = s\mathbf{e}_\phi$ [@malkushydromagnetic1967]. 

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
