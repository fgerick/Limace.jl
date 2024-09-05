### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ b982245e-0558-11ef-0047-8b3236b83ca0
begin
    import Pkg
    Pkg.activate(@__DIR__)
    Pkg.instantiate()
    using CairoMakie, SparseArrays, PlutoUI, GeoMakie
	using Limace.Discretization: discretization_map
end

# ╔═╡ 31cd9b84-eab7-40a4-9150-9f846d25b042
using LinearAlgebra, Limace

# ╔═╡ babb8f6b-6c94-4ef4-8340-0f50f9da7deb
md"""
# Inviscid inertial modes in the sphere

The evolution equation of the velocity ``\mathbf{u}`` is given by the momentum equation

```math
\frac{\partial\mathbf{u}}{\partial t} + 2\Omega \mathbf{e}_z\times\mathbf{u} = -\frac{1}{\rho}\nabla p,
```

satisfying ``\mathbf{u}\cdot\mathbf{n} = 0`` at ``r=1`` (the surface of the sphere).

By projecting this equation onto poloidal and toroidal basis vectors, we eliminate the pressure ``p``.

Assuming ``\mathbf{u}(\mathbf{r},t) = \mathbf{u}(\mathbf{r}) \exp(\lambda t)``, the momentum equation reduces to an eigen problem

```math
\lambda \mathbf{B}\mathbf{x} = \mathbf{A}\mathbf{x},
```

where the eigen vector ``\mathbf{x}`` contain the spectral coefficients of the basis elements and

```math
B_{ij} = \int \mathbf{u}_i \cdot \mathbf{u}_j\,\mathrm{d}V,
```

and

```math
A_{ij} = \int \mathbf{u}_i \cdot \left(2\Omega\mathbf{e}_z\times\mathbf{u}_j\right)\,\mathrm{d}V.
```

For moderate polynomial degrees we can quickly solve this.
"""

# ╔═╡ 2dad3810-855c-47e6-a963-48bc949b75fd
CairoMakie.activate!(; type="png", px_per_unit=1)

# ╔═╡ 1fb533ff-9ed7-4ae7-8a73-1c2fdd8bb161
md"""
## Solve using Limace.jl

First, load the packages:
"""

# ╔═╡ 2edf8beb-5846-4f50-bdb0-79f955d067d9
md"Define a polynomial truncation degree `N`."

# ╔═╡ a863379b-75e5-469d-b580-62ac53fc6426
N = 8

# ╔═╡ 8eb71322-fb38-46c8-bcfe-ca753661a12d
md"Create inviscid velocity basis. Here, we include all azimuthal wave numbers `m = -N:N` (we can also consider each `m` individually for the inviscid inertial modes)."

# ╔═╡ 37fceea1-dbd5-405a-a630-9b99805b0888
basis = Inviscid(N)

# ╔═╡ 7929a5f3-1dec-420a-ab44-6e09a7017e0f
md" Assemble the matrix $\mathbf{A}$ of the Coriolis operator."

# ╔═╡ 1d142fcc-4b20-4026-895c-17cd084b57c3
A = Limace.coriolis(basis)

# ╔═╡ e891d88e-a24b-48c2-97ae-5c448edc59d7
md"The inertial matrix (sometimes called mass matrix) is just the unit operator for the `Inviscid` basis, as it is orthonormal. Due to this orthonormality, the eigen problem reduces to a standard eigenvalue problem. We call `Matrix(A)` to convert the sparse matrix `A` to a dense one. `eigen` gives the dense eigenvalue spectrum. The eigenvalues `λ` have zero real part, i.e. no viscous damping (you can type `λ` by typing `\lambda<tab>`)."

# ╔═╡ 894a4224-3bed-45f5-b884-acd7fe8c16ba
λ, u = eigen(Matrix(A));

# ╔═╡ a9941409-356f-4b30-9cc0-76472d3dad71
md"The eigenvectors are stored in a matrix `u`, so that `λ[i]*u[:,i] ≈ A*u[:,i]`."

# ╔═╡ 2942d573-5c60-41c6-9d72-b79274f322ba
λ[1]*u[:,1] ≈ A*u[:,1]

# ╔═╡ 9c8d542b-40e8-4ffe-ac2f-5803c5eaf59b
md"""
We can compare some of the numerically calculated eigenvalues `λ` to the analytical equation given by Zhang et al. (2008) for equatorially symmetric inertial modes.
"""

# ╔═╡ 34013001-4ff0-4675-bcc3-655e1813bf83
function λ_analytical(m, n) 
	sm = sign(m)
	m = abs(m)
	return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (n * (2n + 2m + 1))) - 1) * im
end;

# ╔═╡ f9e4e190-5b06-40ee-9156-e3d281343e7e
all(m->any(isapprox(λ_analytical(m,1)), λ), 1:N-1)

# ╔═╡ 08f613d1-12d2-4be7-9bb3-a49a80fa6158
md" Since the `Inviscid` basis is complete in the Cartesian polynomials, all values are exactly as the analytical value (up to machine precision). For higher values of `n`, the analytical equation `zhang(m, n)` is only an approximation and one should compare to the roots of the univariate polynomial that gives the eigenvalues of the inertial modes."

# ╔═╡ 878004ba-7462-4595-83d6-b14bb300f653
md"""
## Visualization

We can plot the velocity of the modes, for example on the surface. The azimuthal and latitudinal velocity of the mode corresponding to `zhang(3,1)` is shown below.
"""

# ╔═╡ 8f973c53-9d63-4c6f-b58e-1de7a5ede1fe
nr, nθ, nϕ = 1, 100, 100;

# ╔═╡ c2010494-2d2b-4acd-b09c-a9143cfc470b
r, θ, ϕ = [1.0], (1:nθ)*π/nθ, (0:(nϕ-1))*2π/nϕ;

# ╔═╡ ecde6484-0d34-4789-9588-a0f98db9e0f5
Mr, Mθ, Mϕ = discretization_map(basis,r, θ, ϕ );

# ╔═╡ bcb9d165-c160-48bd-b6b9-d4cefb05450d
imode = findmin(x->abs(x-λ_analytical(3,1)), λ)[2];

# ╔═╡ 6cecaf42-a8c7-40d5-94e5-7a34465d2bdf
# @bind imode PlutoUI.Slider(eachindex(λ), show_value=true)

# ╔═╡ a61b77a4-3e0b-4906-b0c5-fb4d69eada08
_uϕ = real.(reshape(Mϕ*u[:,imode],nθ,nϕ));

# ╔═╡ d1617efb-f212-443f-b8e7-cd26d2769741
_uθ = real.(reshape(Mθ*u[:,imode],nθ,nϕ));

# ╔═╡ 35c4fc9f-e3df-4e0b-9a7d-13db65d5a706
begin
	f = Figure(size=100.0.*(6,2), backgroundcolor = :transparent)
	ax = GeoAxis(f[1,1], dest="proj=moll", yticklabelsvisible=false, xticklabelsvisible=false)
end;

# ╔═╡ 07c55e43-1ce7-4549-b0ec-129512da7dd6
lons, lats = rad2deg.(ϕ).-180,rad2deg.(θ).-90;

# ╔═╡ 52523423-eccc-4ca5-9894-0d1e4ab9ec06
let
	empty!(ax)
	umax = maximum(abs.(_uϕ))
	heatmap!(ax,lons, lats,_uϕ'; colorrange=(-umax,umax), colormap=:balance)
	f
end

# ╔═╡ 68ec3127-f028-43c4-96ee-620a93ce941f
let
	empty!(ax)
	umax = maximum(abs.(_uθ))
	heatmap!(ax,lons, lats,_uθ'; colorrange=(-umax,umax), colormap=:balance)
	f
end

# ╔═╡ Cell order:
# ╟─babb8f6b-6c94-4ef4-8340-0f50f9da7deb
# ╟─b982245e-0558-11ef-0047-8b3236b83ca0
# ╟─2dad3810-855c-47e6-a963-48bc949b75fd
# ╟─1fb533ff-9ed7-4ae7-8a73-1c2fdd8bb161
# ╠═31cd9b84-eab7-40a4-9150-9f846d25b042
# ╟─2edf8beb-5846-4f50-bdb0-79f955d067d9
# ╠═a863379b-75e5-469d-b580-62ac53fc6426
# ╟─8eb71322-fb38-46c8-bcfe-ca753661a12d
# ╠═37fceea1-dbd5-405a-a630-9b99805b0888
# ╟─7929a5f3-1dec-420a-ab44-6e09a7017e0f
# ╠═1d142fcc-4b20-4026-895c-17cd084b57c3
# ╟─e891d88e-a24b-48c2-97ae-5c448edc59d7
# ╠═894a4224-3bed-45f5-b884-acd7fe8c16ba
# ╟─a9941409-356f-4b30-9cc0-76472d3dad71
# ╠═2942d573-5c60-41c6-9d72-b79274f322ba
# ╟─9c8d542b-40e8-4ffe-ac2f-5803c5eaf59b
# ╠═34013001-4ff0-4675-bcc3-655e1813bf83
# ╠═f9e4e190-5b06-40ee-9156-e3d281343e7e
# ╟─08f613d1-12d2-4be7-9bb3-a49a80fa6158
# ╟─878004ba-7462-4595-83d6-b14bb300f653
# ╟─52523423-eccc-4ca5-9894-0d1e4ab9ec06
# ╟─68ec3127-f028-43c4-96ee-620a93ce941f
# ╟─8f973c53-9d63-4c6f-b58e-1de7a5ede1fe
# ╟─c2010494-2d2b-4acd-b09c-a9143cfc470b
# ╟─ecde6484-0d34-4789-9588-a0f98db9e0f5
# ╟─bcb9d165-c160-48bd-b6b9-d4cefb05450d
# ╟─6cecaf42-a8c7-40d5-94e5-7a34465d2bdf
# ╟─a61b77a4-3e0b-4906-b0c5-fb4d69eada08
# ╟─d1617efb-f212-443f-b8e7-cd26d2769741
# ╟─35c4fc9f-e3df-4e0b-9a7d-13db65d5a706
# ╟─07c55e43-1ce7-4549-b0ec-129512da7dd6
