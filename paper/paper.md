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
    affiliation: "1, 2, 3"
affiliations:
 - name: Paris Institute of Earth Physics (IPGP)
   index: 1
   ror: 004gzqz66
 - name: National Centre for Space Studies (CNES), France
   index: 2
   ror: 04h1h0y33
 - name: Royal Observatory of Belgium, Belgium
   index: 3
   ror: 00hjks330
date: 10 June 2024
bibliography: paper.bib
---

# Summary

Hydromagnetic modes in spherical domains are relevant to the liquid cores of planets, moons or stars, as well as rotating fluid dynamics experiments.
These modes are solutions to the linearized rotating magnetohydrodynamic equations that govern electrically conducting fluids under rapid rotation. 
`Limace.jl` is a package written in the Julia programming language [@bezansonjulia2017], based on Galerkin projections of the governing equations onto 
trial vectors of the velocity and magnetic field. It aims to facilitate the calculation of modes in flexible setups with a high-level interface,
whilst remaining computationally performant enough to tackle relevant physical parameters.

# Statement of need

The study of hydromagnetic modes is relevant in particular to Earth's liquid core. 
Despite having been theoretically predicted a long time ago [@hidefree1966; @malkushydromagnetic1967; @braginskytorsional1970], 
recent advances in numerical modelling and new observational evidence in geomagnetic data have reignited interest in these modes [@gerickfast2021; @gilletsatellite2022; @trianacore2022; @luowaves2022a; @luowaves2022; @gerickinterannual2024].
It is therefore relevant to the geophysical and astrophysical fluid dynamics community to have access to a code that models these modes.
The only open-source models to compute hydromagnetic modes in relevant parameters for planetary cores is Kore [@trianaviscous2021], 
a spectral code based on ultraspherical polynomials in axisymmetric setups written in Python.
Although Kore is limited to axisymmetric geometries, the implementation requires substantial work due to the complexity of the underlying spectral equations, resulting in a substantial entry barrier for scientists to model these modes.
In addition to Kore, `Limace.jl` tries to lower this entry barrier, by providing a generic open source model with a very simple high-level API and modern online documentation with practical examples.

Despite having a high-level interface, `Limace.jl` can be used to solve complex and geophysically relevant problems.
A unique feature of `Limace.jl` is the support of complex (non-axisymmetric) background magnetic fields and flows over which the modes evolve.
The code has been developed from the beginning to leave assumptions of symmetry up to the user.
The model code base is tested against mode solutions from the scientific literature to ensure its correctness.

# Theoretical background and implementation details

In order to compute modal solutions, we consider the linearized momentum equation of the incompressible fluid and the linearized induction equation [@iversscalar2008; @gerickinterannual2024]
$$
\begin{aligned}
	\lambda\mathbf{u} =& -\left(\boldsymbol{\nabla}\times\mathbf{u}\right)\times\mathbf{U}_0- \left(\boldsymbol{\nabla}\times\mathbf{U}_0\right)\times\mathbf{u} -2\boldsymbol{\Omega}\times\mathbf{u} - \frac{1}{\rho}\nabla p\\
	 &+ \frac{1}{\rho\mu_0}\left(\left(\boldsymbol{\nabla}\times\mathbf{b}\right)\times\mathbf{B}_0+\left(\boldsymbol{\nabla}\times\mathbf{B}_0\right)\times\mathbf{b}\right) + \nu \boldsymbol{\nabla}^2\mathbf{u},\nonumber\\
	\lambda\mathbf{b} =& \boldsymbol{\nabla}\times\left(\mathbf{U}_0\times\mathbf{b}\right) + \boldsymbol{\nabla}\times\left(\mathbf{u}\times\mathbf{B}_0\right) + \eta \boldsymbol{\nabla}^2\mathbf{b}.
\end{aligned}
$$
with $\mathbf{u}$ the velocity perturbation, $\mathbf{U}_0$ the steady background velocity, $\mathbf{b}$ the magnetic field perturbation, $\mathbf{B}_0$ the background magnetic field, $\boldsymbol{\Omega}$ the rotation axis, $\rho$ the fluid density, $p$ the reduced hydrodynamic pressure, $\mu_0$ the magnetic permeability of free space, $\nu$ the kinematic viscosity, $\eta$ the magnetic diffusivitiy, and $\lambda=-\sigma+\mathrm{i}\omega$, with $\sigma$ the damping rate and $\omega$ the frequency of the oscillatory perturbation to the steady background.

To discretize the linearized equations, they are projected onto trial vectors (Galerkin vector bases) $\mathbf{u}_i$ and $\mathbf{b}_i$ for the velocity and magnetic field, respectively.
Due to the divergence free condition on the velocity and magnetic field, i.e. the flow is incompressible and no magnetic monopoles exist, 
it is convenient to decompose the fields into poloidal and toroidal components, so that
$$
\begin{aligned}
    \mathbf{u} &= \sum_i \alpha_i\mathbf{u}_i = \sum_{l,m,n} \alpha^P_{lmn}\mathbf{P}_{lmn} + \sum_{l,m,n} \alpha^Q_{lmn}\mathbf{Q}_{lmn},\\
    \mathbf{b} &= \sum_i \beta_i\mathbf{b}_i =  \sum_{l,m,n} \beta^S_{lmn}\mathbf{S}_{lmn} + \sum_{l,m,n} \beta^T_{lmn}\mathbf{T}_{lmn},
\end{aligned}
$$
with $\alpha_i, \beta_i \in \mathbb{C}$. 
The respective poloidal $\mathbf{P}_{lmn}$ and toroidal $\mathbf{Q}_{lmn}$ basis vectors in spherical coordinates $(r,\theta,\phi)$ are
$$
\begin{aligned}
	\mathbf{P}_{lmn} & = \boldsymbol{\nabla}\times\boldsymbol{\nabla}\times P_{ln}(r)Y_l^m(\theta,\phi)\mathbf{r}, \\
	\mathbf{Q}_{lmn} & = \boldsymbol{\nabla}\times Q_{ln}(r) Y_l^m(\theta,\phi)\mathbf{r},
\end{aligned}
$$
and analogously for $\mathbf{S}_{lmn}$ and $\mathbf{T}_{lmn}$.
Here, $Y_l^m(\theta,\phi)$ is the (fully normalized) spherical harmonic of degree $l$ and order $m$.
The boundary conditions (or regularity condition at $r=0$) are imposed on the scalar functions $P,S,Q,T$.
The scalar functions can be chosen to have optimal properties, i.e. the resulting basis is orthogonal w.r.t a given inner product [@livermoregalerkin2010; @chenoptimal2018; @gerickinterannual2024].
`Limace.jl` provides several optimal bases that satisfy relevant boundary conditions,
such as non-penetration and no-slip for the velocity and perfectly conducting and insulating exteriors for the magnetic field.

We need to consider all combinations of poloidal and toroidal vector combinations in the projection of the forces.
This leads to several coupling terms, especially for the Lorentz force, advection and induction terms. 
The integrals of these coupling terms over the spherical surfaces are computed through the Adam-Gaunt and Elsasser variables [@jamesadams1973], which are calculated from Wigner symbols (available in Julia through [WignerSymbols.jl](https://github.com/Jutho/WignerSymbols.jl), based on @johanssonfast2016).
The remaining integration in radial direction is done using Gauss-Legendre quadratures, available through [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl).
The exact modelled equations are outlined in @gerickinterannual2024, based on the work of @iversscalar2008.  

From the projected equations, the problem reduces to a generalized eigenproblem
$$
\lambda \mathbf{A}\mathbf{x} = \mathbf{B}\mathbf{x},
$$
that is solved numerically. Here, the eigenvector $\mathbf{x}$ contains the coefficients $\alpha_i$ (and $\beta_i$). The matrix $\mathbf{B}$ is generally not symmetric/Hermitian, but $\mathbf{A}$ can be the unit matrix, symmetric tridiagonal or symmetric, depending on the chosen bases.

For small problem sizes, the eigenproblem can be solved using dense methods, e.g. using the standard library function `LinearAlgebra.eigen`.
For larger problem sizes it is only feasible to compute a few eigensolutions of the sparse system. 
To this end, a shift-invert spectral transform method is provided, based on the sparse LU factorization from `UMFPACK` [@davisalgorithm2004] and the partial Schur decomposition implemented in [ArnoldiMethod.jl](https://github.com/JuliaLinearAlgebra/ArnoldiMethod.jl) [@StoppelsArnoldiMethod].

For postprocessing, `Limace.jl` uses a fast spherical harmonic transform implemented in the [SHTns](https://bitbucket.org/nschaeff/shtns) library [@schaefferefficient2013], and available in Julia through [SHTns.jl](https://github.com/fgerick/SHTns.jl).
It is used to transform the spectral coefficients to vector fields evaluated on a spatial grid.
`Limace.jl` does not provide any plotting routines, but some examples are given, leaving the choice of plotting library up to the user.

# Acknowledgements

I have received funding from the European Research Council (ERC) GRACEFUL Synergy Grant No. 855677. 
This project has been funded by ESA in the framework of EO Science for Society, through contract 4000127193/19/NL/IA (SWARM + 4D Deep Earth: Core). 

I thank Phil Livermore for the key contributions in the theoretical development of the model.

# References
