# Theoretical background

`Limace.jl` uses a Galerkin projection method to discretize a linearized version of the non-linear rotating magnetohydrodynamic (MHD) equations, following works of [bullard_homogeneous_1954, ivers_scalar_2008](@citet) and [gerick_interannual_2024](@citet).
For a very detailed account of all steps, the interested reader is referred to [gerick_interannual_2024](@citet). Here, a quick and simplified overview of the mathematical and physical tools used in the model is given.

In non-linear form, the momentum equations of the incompressible fluid and the induction equation, are written as
```math
\begin{align}
	\frac{\partial \mathbf{U}}{\partial t} + \left(\boldsymbol{\nabla}\times\mathbf{U}\right)\times\mathbf{U} + 2\boldsymbol{\Omega}\times\mathbf{U} &= -\frac{1}{\rho}\nabla P + \frac{1}{\rho\mu_0} \left(\boldsymbol{\nabla}\times\mathbf{B}\right)\times\mathbf{B} + \nu\boldsymbol{\nabla}^2\mathbf{U} + \mathbf{F}, \\
	\frac{\partial \mathbf{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\mathbf{U}\times\mathbf{B}\right) + \eta \boldsymbol{\nabla}^2\mathbf{B},
\end{align}
```
with ``\mathbf{U}`` the velocity, ``\mathbf{B}`` the magnetic field, ``\boldsymbol{\Omega}`` the rotation axis, ``\rho`` the fluid density, ``P`` the reduced hydrodynamic pressure, ``\mu_0`` the magnetic permeability of free space, ``\nu`` the kinematic viscosity, ``\mathbf{F}`` some additional body force, and ``\eta`` the magnetic diffusivitiy.

The velocity, magnetic and pressure fields a linearized, so that 
```math
\begin{align}
	\mathbf{U}(\mathbf{r},t) & = \mathbf{U}_0(\mathbf{r})+ \mathbf{u}(\mathbf{r}) e^{\lambda t}, \\
	\mathbf{B}(\mathbf{r},t) & = \mathbf{B}_0(\mathbf{r})+ \mathbf{b}(\mathbf{r}) e^{\lambda t}, \\
	P(\mathbf{r},t)   & = P_0(\mathbf{r})+ p(\mathbf{r}) e^{\lambda t}.
\end{align}
```
with ``\lambda=-\sigma+\mathrm{i}\omega``, with ``\sigma`` the damping rate and ``\omega`` the frequency of the oscillatory perturbation to the steady background.
Removing the steady part and neglecting higher order terms, the linearized MHD equations read
```math
\begin{align}
	\lambda\mathbf{u} =& -\left(\boldsymbol{\nabla}\times\mathbf{u}\right)\times\mathbf{U}_0- \left(\boldsymbol{\nabla}\times\mathbf{U}_0\right)\times\mathbf{u} -2\Omega\mathbf{e}_z\times\mathbf{u} - \frac{1}{\rho}\nabla p\\
	 &+ \frac{1}{\rho\mu_0}\left(\left(\boldsymbol{\nabla}\times\mathbf{b}\right)\times\mathbf{B}_0+\left(\boldsymbol{\nabla}\times\mathbf{B}_0\right)\times\mathbf{b}\right) + \nu \boldsymbol{\nabla}^2\mathbf{u},\nonumber\\
	\lambda\mathbf{b} =& \boldsymbol{\nabla}\times\left(\mathbf{U}_0\times\mathbf{b}\right) + \boldsymbol{\nabla}\times\left(\mathbf{u}\times\mathbf{B}_0\right) + \eta \boldsymbol{\nabla}^2\mathbf{b}.
\end{align}
```

## Galerkin projection

The linearized equations are projected onto trial vectors ``\boldsymbol{\xi}_i``, so that

```math
f_{ij} = \int \boldsymbol{\xi}_i \cdot \mathbf{f}\left(\mathbf{u}_j,\mathbf{b}_j, \mathbf{U}_0, \mathbf{B}_0\right)\,\mathrm{d}V,
```
where ``\mathbf{f}`` is any of the terms in the considered MHD equations and ``\boldsymbol{\xi}_i = [\mathbf{u}_i, \mathbf{b}_i]``.

## Poloidal and toroidal basis vectors

Due to the divergence free condition on the velocity and magnetic field, i.e. the flow is incompressible and no magnetic monopoles, 
it is convenient to decompose the fields into poloidal and toroidal components.
```math
\begin{align}
    \mathbf{u} &= \sum_i \alpha_i\mathbf{u}_i = \sum_{l,m,n} \alpha^P_{lmn}\mathbf{P}_{lmn} + \sum_{l,m,n} \alpha^Q_{lmn}\mathbf{Q}_{lmn},\\
    \mathbf{b} &= \sum_i \beta_i\mathbf{b}_i =  \sum_{l,m,n} \beta^S_{lmn}\mathbf{S}_{lmn} + \sum_{l,m,n} \beta^T_{lmn}\mathbf{T}_{lmn},
\end{align}
```
with the respective poloidal and toroidal basis vectors
```math
\begin{align}
	\left[\mathbf{P},\mathbf{S}\right]_{lmn} & = \boldsymbol{\nabla}\times\boldsymbol{\nabla}\times \left[P,S\right]_{ln}(r)Y_l^m(\theta,\phi)\mathbf{r}, \\
	\left[\mathbf{Q},\mathbf{T}\right]_{lmn} & = \boldsymbol{\nabla}\times \left[Q,T\right]_{ln}(r) Y_l^m(\theta,\phi)\mathbf{r}.
\end{align}
```
Here, ``Y_l^m`` is the (fully normalized) spherical harmonic of degree ``l`` and order ``m``.

We therefore need to consider all combinations of poloidal and toroidal vector combinations in the projection of the forces.
This leads to several long coupling terms, especially for the Lorentz force and induction term. 
The integrals of these coupling terms over the spherical sourfaces can be reduced to Adam-Gaunt and Elsasser integrals, and it remains to calculate the radial integral. The detailed equations implemented in `Limace.jl` are given in [gerick_interannual_2024](@citet).