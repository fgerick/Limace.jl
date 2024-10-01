# # Inviscid inertial modes in the sphere
#
# The evolution equation of the velocity ``\mathbf{u}`` is given by the momentum equation
#
# ```math
# \frac{\partial\mathbf{u}}{\partial t} + 2\Omega \mathbf{e}_z\times\mathbf{u} = -\frac{1}{\rho}\nabla p,
# ```
#
# satisfying ``\mathbf{u}\cdot\mathbf{n} = 0`` at ``r=1`` (the surface of the sphere).
#
# By projecting this equation onto poloidal and toroidal basis vectors, we eliminate the pressure ``p``.
#
# Assuming ``\mathbf{u}(\mathbf{r},t) = \mathbf{u}(\mathbf{r}) \exp(\lambda t)``, the momentum equation reduces to an eigen problem
#
# ```math
# \lambda \mathbf{B}\mathbf{x} = \mathbf{A}\mathbf{x},
# ```
#
# where the eigen vector ``\mathbf{x}`` contain the spectral coefficients of the basis elements and
#
# ```math
# B_{ij} = \int \mathbf{u}_i \cdot \mathbf{u}_j\,\mathrm{d}V,
# ```
#
# and
#
# ```math
# A_{ij} = \int \mathbf{u}_i \cdot \left(2\Omega\mathbf{e}_z\times\mathbf{u}_j\right)\,\mathrm{d}V.
# ```
#
# For moderate polynomial degrees we can quickly solve this.
#


# ## Solve using Limace.jl
# 
# First, load the packages:

using Limace
using LinearAlgebra, SparseArrays, CairoMakie, GeoMakie
using Limace.Discretization: spectospat
CairoMakie.activate!(; type="png", px_per_unit=2)

# Define a polynomial truncation degree `N`.

N = 8

# Create inviscid velocity basis. Here, we include all azimuthal wave numbers `m = -N:N` (we can also consider each `m` individually for the inviscid inertial modes).

u = Inviscid(N)

# Assemble the matrix $\mathbf{A}$ of the Coriolis operator.

A = Limace.coriolis(u)

# The inertial matrix (sometimes called mass matrix) is just the unit operator for the `Inviscid` basis, as it is orthonormal. Due to this orthonormality, the eigen problem reduces to a standard eigenvalue problem. We call `Matrix(A)` to convert the sparse matrix `A` to a dense one. `eigen` gives the dense eigenvalue spectrum. 
# The eigenvalues `λ` have zero real part, i.e. no viscous damping (you can type `λ` by typing `\lambda<tab>`).

λ, x = eigen(Matrix(A));

# The eigenvectors are stored in a matrix `x`, so that `λ[i]*x[:,i] ≈ A*x[:,i]`.

λ[1]*x[:,1] ≈ A*x[:,1]

# We can compare some of the numerically calculated eigenvalues `λ` to the analytical equation 
# given by [zhang_inertial_2001](@citet) for equatorially symmetric inertial modes.

function λ_analytical(m, n) 
	sm = sign(m)
	m = abs(m)
	return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (n * (2n + 2m + 1))) - 1) * im
end;

all(m->any(isapprox(λ_analytical(m,1)), λ), 1:N-1)

# Since the `Inviscid` basis is complete in the Cartesian polynomials, all values are exactly as the analytical value (up to machine precision). 
# For higher values of `n`, the analytical equation `zhang(m, n)` is only an approximation and one should compare to the roots of the univariate polynomial that gives the eigenvalues of the inertial modes.

#  ## Visualization
#
# We can plot the velocity of the modes, for example on the surface. The azimuthal and latitudinal velocity of the mode corresponding to `zhang(3,1)` is shown below.

r, nθ, nϕ = 1.0, 100, 200

imode = findmin(x->abs(x-λ_analytical(3,1)), λ)[2];

ur, uθ, uϕ, θ, ϕ = spectospat(x[:,imode], u, r, nθ, nϕ);

let
	lons, lats = rad2deg.(ϕ).-180, rad2deg.(θ);
	
	fig = Figure(backgroundcolor = :transparent)
	ax1 = GeoAxis(fig[1,1]; dest="proj=moll", yticklabelsvisible=false, xticklabelsvisible=false, title=L"u_\theta")
	ax2 = GeoAxis(fig[1,2]; dest="proj=moll", yticklabelsvisible=false, xticklabelsvisible=false, title=L"u_\phi")
	_uθ = real.(uθ[1,:,:].+eps())
	umax = maximum(abs,_uθ)
	heatmap!(ax1,lons, lats,_uθ'; colorrange=(-umax,umax), colormap=:balance)
	
	_uϕ = real.(uϕ[1,:,:].+eps())
	umax = maximum(abs,_uϕ)
	heatmap!(ax2,lons, lats,_uϕ'; colorrange=(-umax,umax), colormap=:balance)
	fig
end