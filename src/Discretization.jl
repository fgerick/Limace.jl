module Discretization

using ..Bases
using ..Poly
using ..Poly: dylmdθ, dylmdϕ
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

function discretization_map(u::Basis, r, θ, ϕ)
	nr = length(r)
	nθ = length(θ)
	nϕ = length(ϕ)

	nu = length(u)

	lmnp_u = lmn_p(u)
	lmnt_u = lmn_t(u)

	isr, jsr, aijsr = Int[], Int[], ComplexF64[]
	isθ, jsθ, aijsθ = Int[], Int[], ComplexF64[]
	isϕ, jsϕ, aijsϕ = Int[], Int[], ComplexF64[]


	i = 1
	for r in r, θ in θ, ϕ in ϕ
		j = 1
		for (l,m,n) in lmnp_u
			ur,uθ,uϕ = poloidal_discretize(u, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, ur)
			appendit!(isθ, jsθ, aijsθ, i, j, uθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, uϕ)
			j+=1
		end
		for (l,m,n) in lmnt_u
			ur,uθ,uϕ = toroidal_discretize(u, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, ur)
			appendit!(isθ, jsθ, aijsθ, i, j, uθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, uϕ)
			j+=1
		end
		i+=1
	end


	Mr = sparse(isr, jsr, aijsr, nr*nθ*nϕ,nu+nb)
	Mθ = sparse(isθ, jsθ, aijsθ, nr*nθ*nϕ,nu+nb)
	Mϕ = sparse(isϕ, jsϕ, aijsϕ, nr*nθ*nϕ,nu+nb)
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

	isr, jsr, aijsr = Int[], Int[], ComplexF64[]
	isθ, jsθ, aijsθ = Int[], Int[], ComplexF64[]
	isϕ, jsϕ, aijsϕ = Int[], Int[], ComplexF64[]


	urs = zeros(ComplexF64, )
	i = 1
	for r in r, θ in θ, ϕ in ϕ
		j = 1
		for (l,m,n) in lmnp_u
			ur,uθ,uϕ = poloidal_discretize(u, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, ur)
			appendit!(isθ, jsθ, aijsθ, i, j, uθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, uϕ)
			j+=1
		end
		for (l,m,n) in lmnt_u
			ur,uθ,uϕ = toroidal_discretize(u, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, ur)
			appendit!(isθ, jsθ, aijsθ, i, j, uθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, uϕ)
			j+=1
		end
		for (l,m,n) in lmnp_b
			br,bθ,bϕ = poloidal_discretize(b, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, br)
			appendit!(isθ, jsθ, aijsθ, i, j, bθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, bϕ)
			j+=1
		end
		for (l,m,n) in lmnt_b
			br,bθ,bϕ = toroidal_discretize(b, l,m,n,r,θ,ϕ)
			appendit!(isr, jsr, aijsr, i, j, br)
			appendit!(isθ, jsθ, aijsθ, i, j, bθ)
			appendit!(isϕ, jsϕ, aijsϕ, i, j, bϕ)
			j+=1
		end
		i+=1
	end


	Mr = sparse(isr, jsr, aijsr, nr*nθ*nϕ,nu+nb)
	Mθ = sparse(isθ, jsθ, aijsθ, nr*nθ*nϕ,nu+nb)
	Mϕ = sparse(isϕ, jsϕ, aijsϕ, nr*nθ*nϕ,nu+nb)
	return Mr, Mθ, Mϕ
end

end #module