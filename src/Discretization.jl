module Discretization

using Limace: s,t, lmn_p, lmn_t
using ..Bases
using ..Poly: dylmdθ, dylmdϕ, ∂, ylm
using ..Utils

function poloidal_discretize(::Type{Basis{T}},l,m,n,r,θ,ϕ) where T
	dsrdr = ∂(r->s(Basis{T},l,m,n,r)*r,r)
    ur = l*(l+1)*s(Basis{T}, l,m,n,r)*ylm(l,m,θ,ϕ)/r
    uθ = 1/r*dsrdr*dylmdθ(l,m,θ,ϕ)
    uϕ = 1/(r*sin(θ))*dsrdr*dylmdϕ(l,m,θ,ϕ)
    return (ur,uθ,uϕ)
end

function toroidal_discretize(::Type{Basis{T}},l,m,n,r,θ,ϕ) where T
    ur = 0.0
    uθ = 1/sin(θ)*t(Basis{T}, l,m,n,r)*dylmdϕ(l,m,θ,ϕ)
    uϕ = -t(Basis{T},l,m,n,r)*dylmdθ(l,m,θ,ϕ)

    return (ur,uθ,uϕ)
end

function discretization_map(u::T, r, θ, ϕ) where T<:Basis
	nr = length(r)
	nθ = length(θ)
	nϕ = length(ϕ)

	nu = length(u)

	lmnp_u = lmn_p(u)
	lmnt_u = lmn_t(u)

	Mr = zeros(ComplexF64,nr*nθ*nϕ,nu)
	Mθ = zeros(ComplexF64, nr*nθ*nϕ,nu)
	Mϕ = zeros(ComplexF64, nr*nθ*nϕ,nu)


	i = 1
	for ϕ in ϕ, θ in θ, r in r
		j = 1
		for (l,m,n) in lmnp_u
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = poloidal_discretize(T, l,m,n,r,θ,ϕ)
			j+=1
		end
		for (l,m,n) in lmnt_u
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = toroidal_discretize(T, l,m,n,r,θ,ϕ)
			j+=1
		end
		i+=1
	end


	return Mr, Mθ, Mϕ
end

function discretization_map(u::Basis, b::Basis, r, θ, ϕ)
	nr = length(r)
	nθ = length(θ)
	nϕ = length(ϕ)

	nu = length(u)
	nb = length(b)

	lmnp_u = lmn_p(u)
	lmnt_u = lmn_t(u)

	lmnp_b = lmn_p(b)
	lmnt_b = lmn_t(b)


	Mr = zeros(ComplexF64,nr*nθ*nϕ,nu+nb)
	Mθ = zeros(ComplexF64, nr*nθ*nϕ,nu+nb)
	Mϕ = zeros(ComplexF64, nr*nθ*nϕ,nu+nb)


	i = 1
	for r in r, θ in θ, ϕ in ϕ
		j = 1
		for (l,m,n) in lmnp_u
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = poloidal_discretize(u, l,m,n,r,θ,ϕ)
			j+=1
		end
		for (l,m,n) in lmnt_u
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = toroidal_discretize(u, l,m,n,r,θ,ϕ)
			j+=1
		end
		for (l,m,n) in lmnp_b
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = poloidal_discretize(b, l,m,n,r,θ,ϕ)
			j+=1
		end
		for (l,m,n) in lmnt_b
			Mr[i,j],Mθ[i,j],Mϕ[i,j] = toroidal_discretize(b, l,m,n,r,θ,ϕ)
			j+=1
		end
		i+=1
	end

	return Mr, Mθ, Mϕ
end

end #module