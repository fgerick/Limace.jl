# # Torsional Alfvén modes
#
# Torsional Alfvén modes exist when a conducting fluid is under rapid rotation, so that the rotation time $t_\Omega$ is much shorter than the Alfvén time $t_A=L\sqrt{\mu_0\rho}/B_0$. Put differently, the Lehnert number $\mathrm{Le}= t_\Omega/t_A \ll 1$.
# Then, the flow is dominated by a geostrophic balance at the leading order, and slight departures from the geostrophic balance are restored through the Lorentz force. The result are torsional Alfvén modes (or waves), that are in essence differentially rotating geostrophic cylinders, that stretch the cylindrical radial magnetic field lines that permeate them.
#
# A one-dimensional equation can be derived in the diffusionless limit [braginsky_torsional_1970](@citep).
# Using Limace.jl, we can instead compute weakly diffusive torsional Alfvén modes as solutions to a two-dimensional problem, when the background magnetic field is axisymmetric. In this case, all azimuthal wavenumbers $m$ decouple and we can focus on $m=0$.
#
# In the following example we will solve and illustrate the gravest torsional mode, for an axisymmetric poloidal background magnetic field.


using Limace, SparseArrays, CairoMakie
using Limace.EigenProblem: eigstarget
using Limace.Discretization: spectospat

# ## Assembly

# We can assemble the submatrices and concatenate them to get our final left hand side matrix `LHS` and right hand side matrix `RHS`.
#
# Here, the non-dimensional parameters are the Lehnert number `Le` and the Lundquist number `Lu`. We give the resolution by an integer `N`, that determines the polynomial degree of our solutions.
# This function takes for `B0` a collection of `BasisElement`'s, as shown below.

function assemble(N, Le, Lu, B0)
	u = Inviscid(N; m=0)
	b = Insulating(N; m=0)
	LHS = blockdiag(sparse(Limace.inertial(u),length(u),length(u)), sparse(Limace.inertial(b)))
	RHSc = Limace.coriolis(u)/Le
	RHSl = sum(Limace.lorentz_threaded(u,b,B) for B in B0)
	RHSi = sum(Limace.induction_threaded(b,u,B) for B in B0)
	RHSd = Limace.diffusion(b)/Lu
	RHS = [RHSc RHSl
		   RHSi RHSd]
	return LHS, RHS, u, b
end

# Define parameters and background magnetic field:

N = 100
Le = 5e-4
Lu = 1/Le
B0 = [BasisElement(Basis{Insulating}, Poloidal, (1,0,1),0.3), BasisElement(Basis{Insulating}, Poloidal, (2,0,1),0.7)]

# This background field is a mix of two poloidal field components, `l,m,n=(1,0,1)` and `l,m,n = (2,0,1)`, i.e. dipolar and quadrupolar symmetry.
# 
# Assemble the matrices:

LHS, RHS, u, b = assemble(N, Le, Lu, B0);

# ## Solve the eigen problem
#
# We can calculate some eigenvalues and vectors near a given `target`. For torsional modes it is sensible to choose a target near `λ = -σ+im*ω = 1.0im`, i.e. an Alfvén frequency of 1 (note: this only the case, when using the Alfvén time as the characteristic time-scale in the nondimensionalization).

target = 1.0im

# `Limace.jl` provides the `eigstarget` function, based on shift-invert spectral transform and an implicitly restarted Arnoldi method. The keyword `nev` can be set to the desired number of eigenvalues. 

λ, x = eigstarget(RHS, LHS, target; nev=5);

# ## Plot the solution
#
# We find the solution of largest real part, corresponding to the smallest damping rate. In this simple example this is enough to isolate the gravest torsional mode from the other calculated solutions.

per = sortperm(λ, by = real, rev=true);

λ[per]

# We discretize the corresponding eigenvector in the meridional plane. Since $m=0$, we do not need more than one value of the longitude.

nr, nθ, nϕ = 2N+1, 2N+1, 1

ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(x[:,per[1]], u,b, nr,nθ,nϕ);

# The meridional slices are plotted with `CairoMakie.jl` for example:

let
	f = Figure()
	for (i,(comp, title)) in enumerate(zip((ur, uθ, uϕ),(L"u_r", L"u_\theta", L"u_\phi")))
		ax = Axis(f[1,2i]; aspect=DataAspect(), title)
		hidedecorations!(ax)
		hidespines!(ax)
		x = r .* cos.(θ')
		y = r .* sin.(θ')
		z = real.(comp[:,:,1])
		clim = maximum(abs,z)
		cax = surface!(ax, x,y, z, colormap=:vik, colorrange=(-clim,clim),shading=NoShading)
		Colorbar(f[1,2i-1],cax, flipaxis=false)
	end
	for (i,(comp, title)) in enumerate(zip((br,bθ,bϕ),(L"b_r", L"b_\theta", L"b_\phi")))
		ax = Axis(f[2,2i]; aspect=DataAspect(), title)
		hidedecorations!(ax)
		hidespines!(ax)
		x = r .* cos.(θ')
		y = r .* sin.(θ')
		z = real.(comp[:,:,1])
		clim = maximum(abs,z)
		cax = surface!(ax, x,y, z, colormap=:broc, colorrange=(-clim,clim),shading=NoShading)
		Colorbar(f[2,2i-1],cax, flipaxis=false)
	end

	f
end