### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ 7ce1ff32-7f1d-11ef-1da0-9f9ee589c0f2
let
	import Pkg; Pkg.activate(@__DIR__)
	using Limace, SparseArrays, CairoMakie
	using Limace.EigenProblem: eigstarget
	using Limace.Discretization: spectospat
end

# ╔═╡ d41796ff-8687-405d-b3ad-d97ecddd9f2f
md""" 
# Torsional Alfvén modes

Torsional Alfvén modes exist when a conducting fluid is under rapid rotation, so that the rotation time $t_\Omega$ is much shorter than the Alfvén time $t_A=L\sqrt{\mu_0\rho}/B_0$. Put differently, the Lehnert number $\mathrm{Le}= t_\Omega/t_A \ll 1$.
Then, the flow is dominated by a geostrophic balance at the leading order, and slight departures from the geostrophic balance are restored through the Lorentz force. The result are torsional Alfvén modes (or waves), that are in essence differentially rotating geostrophic cylinders, that stretch the cylindrical radial magnetic field lines that permeate them.

A one-dimensional equation can be derived in the diffusionless limit.
Using Limace.jl, we can instead compute weakly diffusive torsional Alfvén modes as solutions to a two-dimensional problem, when the background magnetic field is axisymmetric. In this case, all azimuthal wavenumbers $m$ decouple and we can focus on $m=0$.

In the following example we will solve and illustrate the gravest torsional mode, for an axisymmetric poloidal background magnetic field.

"""

# ╔═╡ 40b49be7-8d6b-4b5f-91eb-a7e759ce52a1
md"""
## Assembly

We can assemble the submatrices and concatenate them to get our final left hand side matrix `LHS` and right hand side matrix `RHS`.

Here, the non-dimensional parameters are the Lehnert number `Le` and the Lundquist number `Lu`. We give the resolution by an integer `N`, that determines the polynomial degree of our solutions.
This function takes for `B0` a collection of `BasisElement`'s, as shown below.
"""

# ╔═╡ 9624ec91-79a5-4730-bba2-66c676570fb6
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

# ╔═╡ 40077073-8d10-4c5a-977e-4d91df3b234e
md"""
Define parameters and background magnetic field:
"""

# ╔═╡ 2b989642-6eb5-4396-b72a-3fd655124386
begin
	N = 100
	Le = 5e-4
	Lu = 1/Le
	B0 = [BasisElement(Basis{Insulating}, Poloidal, (1,0,1),0.3), BasisElement(Basis{Insulating}, Poloidal, (2,0,1),0.7)]
end

# ╔═╡ 5083c4cb-4837-4a75-b506-e24b0e41f9e1
md"""

This background field is a mix of two poloidal field components, `l,m,n=(1,0,1)` and `l,m,n = (2,0,1)`, i.e. dipolar and quadrupolar symmetry.

Assemble the matrices:
"""

# ╔═╡ 5878483c-56c6-458a-b577-6926e8f1942a
LHS, RHS, u, b = assemble(N, Le, Lu, B0);

# ╔═╡ 29882d33-5f8c-411b-a8c0-d94d14db0f4a
RHS

# ╔═╡ 8dcd39ec-2a5a-49e3-8dcb-bc04d1464a2d
LHS

# ╔═╡ 07c9fe67-fb95-40cf-a9fb-7d4faa780c63
md"""
## Solve the eigen problem

We can calculate some eigenvalues and vectors near a given `target`. For torsional modes it is sensible to choose a target near `λ = -σ+im*ω = 1.0im`, i.e. an Alfvén frequency of 1 (note: this only the case, when using the Alfvén time as the characteristic time-scale in the nondimensionalization).

`Limace.jl` provides the `eigstarget` function, based on shift-invert spectral transform and an implicitly restarted Arnoldi method. The keyword `nev` can be set to the desired number of eigenvalues. 
"""

# ╔═╡ 600fd91e-9fec-463b-8c1a-6d6611fb55c2
target = 1.0im

# ╔═╡ 6e815c61-386d-415a-86e3-4219cce99685
λ, x = eigstarget(RHS, LHS, target; nev=5);

# ╔═╡ bc3f1b76-3ef7-4fef-a6f6-b201b8e67140
md"""
## Plot the solution

We find the solution of largest real part, corresponding to the smallest damping rate. In this simple example this is enough to isolate the gravest torsional mode from the other calculated solutions.
"""

# ╔═╡ 62475ddc-571b-4c68-a0db-54a97fc4a9d4
per = sortperm(λ, by = real, rev=true);

# ╔═╡ 5176e295-0375-4f83-8207-109eb58a0fa1
λ[per]

# ╔═╡ ae61d632-56af-4350-ac57-24440f2bfc89
md"""
We discretize the corresponding eigenvector in the meridional plane. Since $m=0$, we do not need more than one value of the longitude.
"""

# ╔═╡ b6f7b17d-b3ec-4f8b-a792-ef8ec73956e3
nr, nθ, nϕ = 2N+1, 2N+1, 1

# ╔═╡ 04f2a755-31ea-42ab-9103-0a237de21f6e
ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(x[:,per[1]], u,b, nr,nθ,nϕ);

# ╔═╡ 75b400a5-47b7-4529-a926-def061bed055
md"""
The meridional slices are plotted with `CairoMakie.jl` for example:
"""

# ╔═╡ d4db0ab5-b986-4629-9fdc-1932e383404c
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

# ╔═╡ Cell order:
# ╟─d41796ff-8687-405d-b3ad-d97ecddd9f2f
# ╠═7ce1ff32-7f1d-11ef-1da0-9f9ee589c0f2
# ╟─40b49be7-8d6b-4b5f-91eb-a7e759ce52a1
# ╠═9624ec91-79a5-4730-bba2-66c676570fb6
# ╟─40077073-8d10-4c5a-977e-4d91df3b234e
# ╠═2b989642-6eb5-4396-b72a-3fd655124386
# ╟─5083c4cb-4837-4a75-b506-e24b0e41f9e1
# ╠═5878483c-56c6-458a-b577-6926e8f1942a
# ╠═29882d33-5f8c-411b-a8c0-d94d14db0f4a
# ╠═8dcd39ec-2a5a-49e3-8dcb-bc04d1464a2d
# ╠═07c9fe67-fb95-40cf-a9fb-7d4faa780c63
# ╠═600fd91e-9fec-463b-8c1a-6d6611fb55c2
# ╠═6e815c61-386d-415a-86e3-4219cce99685
# ╟─bc3f1b76-3ef7-4fef-a6f6-b201b8e67140
# ╠═62475ddc-571b-4c68-a0db-54a97fc4a9d4
# ╠═5176e295-0375-4f83-8207-109eb58a0fa1
# ╟─ae61d632-56af-4350-ac57-24440f2bfc89
# ╠═b6f7b17d-b3ec-4f8b-a792-ef8ec73956e3
# ╠═04f2a755-31ea-42ab-9103-0a237de21f6e
# ╟─75b400a5-47b7-4529-a926-def061bed055
# ╠═d4db0ab5-b986-4629-9fdc-1932e383404c
