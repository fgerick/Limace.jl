# # Torsional Alfvén modes
#
# Torsional Alfvén modes exist when a conducting fluid is under rapid rotation, so that the rotation time $t_\Omega$ is much shorter than the Alfvén time $t_A=L\sqrt{\mu_0\rho}/B_0$. 
# Put differently, the Lehnert number $\mathrm{Le}= t_\Omega/t_A \ll 1$.
# Then, the flow is dominated by a geostrophic balance at the leading order, and slight departures from the geostrophic balance are restored through the Lorentz force. 
# The result are torsional Alfvén modes (or waves), that are in essence differentially rotating geostrophic cylinders, that stretch the cylindrical radial magnetic field lines that permeate them.
#
# A one-dimensional equation can be derived in the diffusionless limit [braginsky_torsional_1970](@citep).
# Using Limace.jl, we can instead compute weakly diffusive torsional Alfvén modes as solutions to a two-dimensional problem, when the background magnetic field is axisymmetric. 
# In this case, all azimuthal wavenumbers $m$ decouple and we can focus on $m=0$.
#
# In the following example we will solve and illustrate the gravest torsional mode, for an axisymmetric poloidal background magnetic field.


using Limace, LinearAlgebra, SparseArrays, CairoMakie
using Limace.EigenSolve: eigstarget
using Limace.Discretization: spectospat

# ## Assembly

# We can assemble the submatrices and concatenate them to get our final left hand side matrix `LHS` and right hand side matrix `RHS`.
#
# Here, the non-dimensional parameters are the Lehnert number `Le` and the Lundquist number `Lu`. We give the resolution by an integer `N`, that determines the polynomial degree of our solutions.
# We define `B0` as a collection of two `BasisElement`'s, an arbitrary mix of two poloidal field components, `l,m,n=(1,0,1)` and `l,m,n = (2,0,1)`, i.e. dipolar and quadrupolar symmetry.
# Define parameters, inviscid and insulating bases, and background magnetic field:

N = 80
Le = 5e-4
Lu = 1/Le

u = Inviscid(N; m=0)
b = Insulating(N; m=0)
B0 = [BasisElement(b, Poloidal, (1,0,1),0.3), BasisElement(b, Poloidal, (2,0,1),0.7)]

# The factors `0.3` and `0.7` are arbitrary, but they can be used to adjust the relative strength of the components of the background magnetic field.
# We can now put it together in a `LimaceProblem`.
bases = [u, b]
forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u,b,B0), Limace.InductionB0(b,u,B0), Limace.Diffusion(b, 1/Lu)]
problem = LimaceProblem(bases, forcings);

# Assemble the problem matrices:

Limace.assemble!(problem; threads=true);

# ## Solve the eigen problem
#
# We can calculate some eigenvalues and vectors near a given `target`. For torsional modes it is sensible to choose a target near `λ = -σ+im*ω = 1.0im`, 
# i.e. an Alfvén frequency of 1 (note: this only the case, when using the Alfvén time as the characteristic time-scale in the nondimensionalization).

target = 1.0im

# The problem can be solved using `Limace.solve!` with the `:sparse` method and the `target` keyword argument, which is based on shift-invert spectral transform and an implicitly restarted Arnoldi method. 
# The keyword `nev` can be set to the desired number of eigenvalues.

Limace.solve!(problem; method=:sparse, target, nev=5);

# The `LimaceProblem` object now contains the eigenvalues and eigenvectors in `problem.sol.values` and `problem.sol.vectors`, respectively. 

λ, x = problem.sol.values, problem.sol.vectors;

# Alternatively, one can use the `eigstarget` function using the problem matrices directly:

λ, x = eigstarget(problem.RHS, problem.LHS, target; nev=5);

# ## Plot the solution
#
# We find the solution of largest real part, corresponding to the smallest damping rate. 
# In this simple example this is enough to isolate the gravest torsional mode from the other calculated solutions.

per = sortperm(λ, by = real, rev=true);
λ[per]

# We discretize the corresponding eigenvector in the meridional plane. Since $m=0$, we do not need more than one value of the longitude.

nr, nθ, nϕ = 2N+1, 2N+1, 1
y = x[:,per[1]]
y /= y[Limace.Bases.lmn2k_t_dict(u)[(1,0,1)] + Limace.Bases.np(u)] #norm by m=0 toroidal component.
ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(y, u,b, nr,nθ,nϕ);

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