module Discretization

using Limace: s, t, lmn_p, lmn_t
using ..Bases
using ..Bases: Sphere
using ..Poly: dylmdθ, dylmdϕ, ∂, ylm
using ..Quadrature
using StaticArrays
using SHTns


function poloidal_discretize(::Type{Basis{T}}, V::Volume, l, m, n, r, θ, ϕ) where {T}
    if r <= V.r1
        dsrdr = ∂(r -> s(Basis{T}, V, l, m, n, r) * r, r)
        ur = l * (l + 1) * s(Basis{T}, V, l, m, n, r) * ylm(l, m, θ, ϕ) / r
        uθ = 1 / r * dsrdr * dylmdθ(l, m, θ, ϕ)
        uϕ = 1 / (r * sin(θ)) * dsrdr * dylmdϕ(l, m, θ, ϕ)
        return SVector(ur, uθ, uϕ)
    else
        sc = s(Basis{T}, V, l, m, n, V.r1)
        ur = l * sc * ylm(l, m, θ, ϕ) * (l + 1) * r^(-l - 2)
        uθ = -l / r * sc * r^(-l - 1) * dylmdθ(l, m, θ, ϕ)
        uϕ = -l / (r * sin(θ)) * sc * dylmdϕ(l, m, θ, ϕ)
        return SVector(ur, uθ, uϕ)
    end
end


function toroidal_discretize(::Type{Basis{T}}, V::Volume, l, m, n, r, θ, ϕ) where {T}
    if r <= V.r1
        ur = 0.0
        uθ = 1 / sin(θ) * t(Basis{T}, V, l, m, n, r) * dylmdϕ(l, m, θ, ϕ)
        uϕ = -t(Basis{T}, V, l, m, n, r) * dylmdθ(l, m, θ, ϕ)
        return SVector(ur, uθ, uϕ)
    else
        return SVector(0.0, 0.0, 0.0)
    end
end

function discretize(b::BasisElement{TB,TP,T}, r, θ, ϕ, V::Volume=Sphere()) where {TB<:Basis,TP<:Helmholtz,T<:Number}
    if TP <: Poloidal
        return b.factor * poloidal_discretize(TB, V, b.lmn..., r, θ, ϕ)
    else
        return b.factor * toroidal_discretize(TB, V, b.lmn..., r, θ, ϕ)
    end
end

function discretize(bs, r, θ, ϕ, V::Volume=Sphere())
    return mapreduce(b -> discretize(b, r, θ, ϕ, V), +, bs)
end

# function discretize(αs::Vector{T}, u::TU, r, θ, ϕ) where {T<:Number,TU<:Basis}
#     @assert length(αs) == length(u)
#     nr = length(r)
#     nθ = length(θ)
#     nϕ = length(ϕ)

#     lmnp_u = lmn_p(u)
#     lmnt_u = lmn_t(u)

#     ur = zeros(ComplexF64, nr, nθ, nϕ)
#     uθ = zeros(ComplexF64, nr, nθ, nϕ)
#     uϕ = zeros(ComplexF64, nr, nθ, nϕ)

#     _np = length(lmnp_u)
#     αsp = @view αs[1:_np]
#     αst = @view αs[_np+1:end]

#     lck = ReentrantLock()

#     Threads.@threads for _i in eachindex(lmnp_u)
#         (l, m, n) , α= lmnp_u[_i], αsp[_i]
#         @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
#             _ur,_uθ,_uϕ = poloidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
#             Threads.lock(lck) do
#                 ur[i, j, k] += α*_ur
#                 uθ[i, j, k] += α*_uθ
#                 uϕ[i, j, k] += α*_uϕ
#             end
#         end
#     end
#     Threads.@threads for _i in eachindex(lmnt_u)
#         (l, m, n) ,α = lmnt_u[_i], αst[_i]
#         @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
#             _ur,_uθ,_uϕ = toroidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
#             Threads.lock(lck) do
#                 ur[i, j, k] += α*_ur
#                 uθ[i, j, k] += α*_uθ
#                 uϕ[i, j, k] += α*_uϕ
#             end
#         end
#     end


#     return ur, uθ, uϕ
# end

function discretize(αs::Vector{T}, u::TU, rs, θ, ϕ) where {T<:Number,TU<:Basis}
    @assert length(αs) == length(u)
    nr = length(rs)
    nθ = length(θ)
    nϕ = length(ϕ)

    lmnp_u = lmn_p(u)
    lmnt_u = lmn_t(u)

    ur = zeros(ComplexF64, nr, nθ, nϕ)
    uθ = zeros(ComplexF64, nr, nθ, nϕ)
    uϕ = zeros(ComplexF64, nr, nθ, nϕ)

    _np = length(lmnp_u)
    αsp = @view αs[1:_np]
    αst = @view αs[_np+1:end]


    Threads.@threads for i in eachindex(rs)
        r = rs[i]
        for ((l, m, n) , α) in zip(lmnp_u,αsp), (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ)
            _ur,_uθ,_uϕ = poloidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
            ur[i, j, k] += α*_ur
            uθ[i, j, k] += α*_uθ
            uϕ[i, j, k] += α*_uϕ
        end
        for ((l,m,n), α) in zip(lmnt_u, αst), (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ)
            _ur,_uθ,_uϕ = toroidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
            ur[i, j, k] += α*_ur
            uθ[i, j, k] += α*_uθ
            uϕ[i, j, k] += α*_uϕ
        end
    end


    return ur, uθ, uϕ
end

function discretize(αs::Vector{T}, u::TU, b::TB, r, θ, ϕ) where {T<:Number,TU<:Basis, TB<:Basis}
    @assert length(αs) == length(u)+length(b)
    nr = length(r)
    nθ = length(θ)
    nϕ = length(ϕ)

    #velocity

    nu = length(u)
    lmnp_u = lmn_p(u)
    lmnt_u = lmn_t(u)

    ur = zeros(ComplexF64, nr, nθ, nϕ)
    uθ = zeros(ComplexF64, nr, nθ, nϕ)
    uϕ = zeros(ComplexF64, nr, nθ, nϕ)

    _np = length(lmnp_u)
    αspu = @view αs[1:_np]
    αstu = @view αs[_np+1:nu]

    lck = ReentrantLock()

    Threads.@threads for _i in eachindex(lmnp_u)
        (l, m, n) , α= lmnp_u[_i], αspu[_i]
        @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
            _ur,_uθ,_uϕ = poloidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
            Threads.lock(lck) do
                ur[i, j, k] += α*_ur
                uθ[i, j, k] += α*_uθ
                uϕ[i, j, k] += α*_uϕ
            end
        end
    end
    Threads.@threads for _i in eachindex(lmnt_u)
        (l, m, n) ,α = lmnt_u[_i], αstu[_i]
        @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
            _ur,_uθ,_uϕ = toroidal_discretize(TU, u.V, l, m, n, r, θ, ϕ)
            Threads.lock(lck) do
                ur[i, j, k] += α*_ur
                uθ[i, j, k] += α*_uθ
                uϕ[i, j, k] += α*_uϕ
            end
        end
    end


    #mag. field
    
    lmnp_b = lmn_p(b)
    lmnt_b = lmn_t(b)

    br = zeros(ComplexF64, nr, nθ, nϕ)
    bθ = zeros(ComplexF64, nr, nθ, nϕ)
    bϕ = zeros(ComplexF64, nr, nθ, nϕ)

    _npb = length(lmnp_b)
    αspb = @view αs[nu+1:nu+_npb]
    αstb = @view αs[nu+_npb+1:end]

    lck = ReentrantLock()

    Threads.@threads for _i in eachindex(lmnp_b)
        (l, m, n) , α= lmnp_b[_i], αspb[_i]
        @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
            _ur,_uθ,_uϕ = poloidal_discretize(TB, b.V, l, m, n, r, θ, ϕ)
            Threads.lock(lck) do
                br[i, j, k] += α*_ur
                bθ[i, j, k] += α*_uθ
                bϕ[i, j, k] += α*_uϕ
            end
        end
    end
    Threads.@threads for _i in eachindex(lmnt_b)
        (l, m, n) ,α = lmnt_b[_i], αstb[_i]
        @inbounds for (k,ϕ) in enumerate(ϕ), (j,θ) in enumerate(θ), (i,r) in enumerate(r)
            _ur,_uθ,_uϕ = toroidal_discretize(TB, b.V, l, m, n, r, θ, ϕ)
            Threads.lock(lck) do
                br[i, j, k] += α*_ur
                bθ[i, j, k] += α*_uθ
                bϕ[i, j, k] += α*_uϕ
            end
        end
    end

    return ur, uθ, uϕ, br, bθ, bϕ
end

function discretization_map(u::T, r, θ, ϕ) where {T<:Basis}
    nr = length(r)
    nθ = length(θ)
    nϕ = length(ϕ)

    nu = length(u)

    lmnp_u = lmn_p(u)
    lmnt_u = lmn_t(u)

    Mr = zeros(ComplexF64, nr * nθ * nϕ, nu)
    Mθ = zeros(ComplexF64, nr * nθ * nϕ, nu)
    Mϕ = zeros(ComplexF64, nr * nθ * nϕ, nu)


    i = 1
    for ϕ in ϕ, θ in θ, r in r
        j = 1
        for (l, m, n) in lmnp_u
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = poloidal_discretize(T, u.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        for (l, m, n) in lmnt_u
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = toroidal_discretize(T, u.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        i += 1
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


    Mr = zeros(ComplexF64, 2nr * nθ * nϕ, nu + nb)
    Mθ = zeros(ComplexF64, 2nr * nθ * nϕ, nu + nb)
    Mϕ = zeros(ComplexF64, 2nr * nθ * nϕ, nu + nb)


    i = 1
    for r in r, θ in θ, ϕ in ϕ
        j = 1
        for (l, m, n) in lmnp_u
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = poloidal_discretize(typeof(u), u.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        for (l, m, n) in lmnt_u
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = toroidal_discretize(typeof(u), u.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        i += 1
    end
    for r in r, θ in θ, ϕ in ϕ
        j = nu + 1
        for (l, m, n) in lmnp_b
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = poloidal_discretize(typeof(b), b.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        for (l, m, n) in lmnt_b
            Mr[i, j], Mθ[i, j], Mϕ[i, j] = toroidal_discretize(typeof(b), b.V, l, m, n, r, θ, ϕ)
            j += 1
        end
        i += 1
    end

    return Mr, Mθ, Mϕ
end


function coeffs_to_SHTnSlmTlm!(Qlmu, Slmu, Tlmu, Qlmb, Slmb, Tlmb, lmnpu, lmntu, lmnpb, lmntb, coeffs, sht, u::TU, b::TB, r) where {TU<:Basis, TB<:Basis}

	npu = length(lmnpu)
	ntu = length(lmntu)
	nu = npu+ntu
	npb = length(lmnpb)
	ntb = length(lmntb)

    # @show length(coeffs)
    # @show nu + npb + ntb
    # @assert length(coeffs) == nu + npb +ntb

	#velocity

	for (c,(l,m,n)) in zip(@views(coeffs[1:npu]),lmnpu)
		@inline p = r-> s(TU,u.V,l,m,n,r)
		Qlmu[SHTns.LM_cplx(sht,l,m)] += c*l*(l+1)/r*p(r)
		Slmu[SHTns.LM_cplx(sht,l,m)] += c/r*∂(r->r*p(r),r)

	end
	for (c,(l,m,n)) in zip(@views(coeffs[npu+1:npu+ntu]),lmntu)
		Tlmu[SHTns.LM_cplx(sht,l,m)] += c*t(TU,u.V,l,m,n,r)
	end

	#magnetic field

	for (c,(l,m,n)) in zip(@views(coeffs[nu+1:nu+npb]),lmnpb)
		@inline p = r-> s(TB,b.V,l,m,n,r)
		Qlmb[SHTns.LM_cplx(sht,l,m)] += c*l*(l+1)/r*p(r)
		Slmb[SHTns.LM_cplx(sht,l,m)] += c/r*∂(@inline(r->r*p(r)),r)

	end
	for (c,(l,m,n)) in zip(@views(coeffs[nu+npb+1:nu+npb+ntb]),lmntb)
		Tlmb[SHTns.LM_cplx(sht,l,m)] += c*t(TB,b.V,l,m,n,r)
	end

	return nothing

end

function coeffs_to_SHTnSlmTlm!(Qlmu, Slmu, Tlmu, lmnpu, lmntu, coeffs, sht, u::TU, r) where {TU<:Basis}

	npu = length(lmnpu)
	ntu = length(lmntu)
	nu = npu+ntu
    @assert length(coeffs) == nu

	#velocity

	for (c,(l,m,n)) in zip(@views(coeffs[1:npu]),lmnpu)
		@inline p = r-> s(TU,u.V,l,m,n,r)
		Qlmu[SHTns.LM_cplx(sht,l,m)] += c*l*(l+1)/r*p(r)
		Slmu[SHTns.LM_cplx(sht,l,m)] += c/r*∂(r->r*p(r),r)

	end
	for (c,(l,m,n)) in zip(@views(coeffs[npu+1:npu+ntu]),lmntu)
		Tlmu[SHTns.LM_cplx(sht,l,m)] += c*t(TU,u.V,l,m,n,r)
	end

	return nothing

end

function discretizationsetup(L::Int, V::Volume)
	
	nr = 2L
	sht = SHTnsCfg(L)
	lats, lons = SHTns.grid(sht; colat=true)
	rgrid, _ = Quadrature.rquad(nr,V.r0, V.r1)
	return sht, rgrid, lats, lons
end

function discretizationsetup(L::Int, V::Volume, nr, nlat, nlon; mmax=L)
	
	sht = SHTnsCfg(L, mmax, 1, nlat, nlon)
	lats, lons = SHTns.grid(sht; colat=true)
	rgrid, _ = Quadrature.rquad(nr,V.r0, V.r1)
	return sht, rgrid, lats, lons
end

function discretizationsetup(L::Int, nlat, nlon; mmax=L)
	sht = SHTnsCfg(L, mmax, 1, nlat, nlon)
	lats, lons = SHTns.grid(sht; colat=true)
	return sht, lats, lons
end

function setup_coeffs_to_SHTnSlmTlm(sht, u::TU, b::TB) where {TU<:Basis, TB<:Basis}
	lmnpu = lmn_p(u)
	lmntu = lmn_t(u)
	lmnpb = lmn_p(b)
	lmntb = lmn_t(b)
	

	#velocity
	Qlmu = zeros(ComplexF64,sht.nlm_cplx)
	Slmu = zeros(ComplexF64,sht.nlm_cplx)
	Tlmu = zeros(ComplexF64,sht.nlm_cplx)

	#magnetic field
	Qlmb = zeros(ComplexF64,sht.nlm_cplx)
	Slmb = zeros(ComplexF64,sht.nlm_cplx)
	Tlmb = zeros(ComplexF64,sht.nlm_cplx)

	return lmnpu, lmntu, lmnpb, lmntb, Qlmu, Slmu, Tlmu, Qlmb, Slmb, Tlmb
	

end

function setup_coeffs_to_SHTnSlmTlm(sht, u::TU) where {TU<:Basis}
	lmnpu = lmn_p(u)
	lmntu = lmn_t(u)
	
	Qlmu = zeros(ComplexF64,sht.nlm_cplx)
	Slmu = zeros(ComplexF64,sht.nlm_cplx)
	Tlmu = zeros(ComplexF64,sht.nlm_cplx)

	return lmnpu, lmntu, Qlmu, Slmu, Tlmu
	

end

function _spectospat_shtns(coeffs::Vector{T}, sht, rgrid, u::TU) where {TU<:Basis, T<:ComplexF64}
	nr = length(rgrid)
	nϕ = sht.nphi
	nθ = sht.nlat
	ur = zeros(T,nr,nθ,nϕ)
	uθ = zeros(T,nr,nθ,nϕ)
	uϕ = zeros(T,nr,nθ,nϕ)
	_ur = zeros(T,nθ,nϕ)
	_uθ, _uϕ = similar(_ur), similar(_ur)
	lmnpu, lmntu, Qlmu, Slmu, Tlmu = setup_coeffs_to_SHTnSlmTlm(sht, u)
	for (ir, r) in enumerate(rgrid)
		coeffs_to_SHTnSlmTlm!(Qlmu, Slmu, Tlmu, lmnpu, lmntu, coeffs, sht, u, r)
		SHTns.synth!(sht, Qlmu, Slmu, Tlmu, _ur, _uθ, _uϕ)
		
		@inbounds for j in axes(_ur,2), i in axes(_ur,1)
			ur[ir,i,j] = _ur[i,j]
			uθ[ir,i,j] = _uθ[i,j]
			uϕ[ir,i,j] = _uϕ[i,j]
		end
		@. Qlmu = Slmu = Tlmu = 0.0
	end
	return ur,uθ,uϕ

end

function _spectospat_shtns(coeffs::Vector{T}, sht, rgrid, u::TU, b::TB) where {TU<:Basis, TB<:Basis, T<:ComplexF64}
	nr = length(rgrid)
	nϕ = sht.nphi
	nθ = sht.nlat
	ur = zeros(T,nr,nθ,nϕ)
	uθ = zeros(T,nr,nθ,nϕ)
	uϕ = zeros(T,nr,nθ,nϕ)
	br = zeros(T,nr,nθ,nϕ)
	bθ = zeros(T,nr,nθ,nϕ)
	bϕ = zeros(T,nr,nθ,nϕ)
	_ur = zeros(T,nθ,nϕ)
	_uθ, _uϕ, _br, _bϕ, _bθ = similar(_ur), similar(_ur), similar(_ur), similar(_ur), similar(_ur)
	lmnpu, lmntu, lmnpb, lmntb, Qlmu, Slmu, Tlmu, Qlmb, Slmb, Tlmb = setup_coeffs_to_SHTnSlmTlm(sht, u, b)
	for (ir, r) in enumerate(rgrid)
		coeffs_to_SHTnSlmTlm!(Qlmu, Slmu, Tlmu, Qlmb, Slmb, Tlmb, lmnpu, lmntu, lmnpb, lmntb, coeffs, sht, u, b, r)
		SHTns.synth!(sht, Qlmu, Slmu, Tlmu, _ur, _uθ, _uϕ)
		SHTns.synth!(sht, Qlmb, Slmb, Tlmb, _br, _bθ, _bϕ)
		
		for j in axes(_ur,2), i in axes(_ur,1)
			ur[ir,i,j] = _ur[i,j]
			uθ[ir,i,j] = _uθ[i,j]
			uϕ[ir,i,j] = _uϕ[i,j]
			br[ir,i,j] = _br[i,j]
			bθ[ir,i,j] = _bθ[i,j]
			bϕ[ir,i,j] = _bϕ[i,j]
		end
		@. Qlmu = Slmu = Tlmu = Qlmb = Slmb = Tlmb = 0.0
	end
	return ur,uθ,uϕ, br,bθ,bϕ

end

function spectospat(coeffs::Vector{T}, u::TU, nr::Int, nθ::Int, nϕ::Int) where {TU<:Basis, T<:ComplexF64}
    sht, r, θ, ϕ = discretizationsetup(u.N, u.V, nr, nθ, nϕ; mmax=maximum(u.m))
    ur,uθ,uϕ = _spectospat_shtns(coeffs, sht, r, u)
    return ur,uθ,uϕ, r,θ,ϕ
end


function spectospat(coeffs::Vector{T}, u::TU, r::Float64, nθ::Int, nϕ::Int) where {TU<:Basis, T<:ComplexF64}
    sht, θ, ϕ = discretizationsetup(u.N, nθ, nϕ; mmax=maximum(u.m))
    ur,uθ,uϕ = _spectospat_shtns(coeffs, sht, r, u)
    return ur,uθ,uϕ, θ,ϕ
end

function spectospat(coeffs::Vector{T}, u::TU, b::TB, nr::Int, nθ::Int, nϕ::Int) where {TU<:Basis, TB<:Basis, T<:ComplexF64}
    sht, r, θ, ϕ = discretizationsetup(u.N, u.V, nr, nθ, nϕ; mmax=maximum(u.m))
    ur,uθ,uϕ, br,bθ,bϕ = _spectospat_shtns(coeffs, sht, r, u, b)
    return ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ
end


function spectospat(coeffs::Vector{T}, u::TU, b::TB, r::Tr, nθ::Int, nϕ::Int) where {TU<:Basis, TB<:Basis, T<:ComplexF64, Tr<:Union{AbstractVector{Float64},Float64}}
    sht, θ, ϕ = discretizationsetup(u.N, nθ, nϕ; mmax=maximum(u.m))
    ur,uθ,uϕ, br,bθ,bϕ = _spectospat_shtns(coeffs, sht, r, u, b)
    return ur,uθ,uϕ, br,bθ,bϕ, θ,ϕ
end


end #module
