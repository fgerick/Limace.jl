
@testset "unit sphere norm of n=1 poloidal insulating MF fields" begin
	N = 20
	r,wr = Limace.Quadrature.rquad(N)
	b = Limace.Insulating(N)

	for l in 1:5, m in -l:l, n in 1:3
		@test Limace._inertial_ss(b, (l,m,n), (l,m,n), r,wr; external=false)*Limace.InsulatingBasis.unitspherenorm(l,n)^2 ≈ 1.0
		# @test Limace._inertial_SS((l,m,n), (l,m,n), r,wr, DP.s_mf,DP.s_mf)*Limace.InsulatingMFBasis.unitspherenorm(l,n)^2 ≈ 1.0
	end
end

@testset "normalizing non-orthogonal background fields so that total (4π/3)⁻¹ ∫₀¹  B⋅B dV = 1" begin

	using Limace.InsulatingBasis: unitspherenorm, norm_B0fac!, inner_b0norm

	nf = unitspherenorm

	Bs = [( [(1,0,1)], [], [nf(1,1)])]
	push!(Bs, ( [(1,0,1)], [], [1.0]))
	push!(Bs, ( [(1,0,1),(1,0,2)],[], [1.0,1.0]))
	push!(Bs, ( [(1,0,1),(1,0,2),(2,0,1)],[], [1.0,1.0,0.5]))
	push!(Bs, ( [(1,0,1),(1,0,2), (2,0,1)],[], [1.0,1.0,1.0]))
	push!(Bs, ([(1,0,1),(1,0,2), (2,0,1)],[],  [1.0,9.0,1.0]))
	push!(Bs, ( [(1,0,1),(1,0,2), (2,0,1)],[(1,0,1)], [1.0,1.0,1.0,1.0]))
	push!(Bs, ( [(1,0,1),(1,0,2), (2,0,1)],[(1,0,1)], [1,1,1,1.0]))
	
	for (lmn_p, lmn_t, B0fac) in Bs
		norm_B0fac!(B0fac, lmn_p, lmn_t )
		A = inner_b0norm(lmn_p, lmn_t)
		@test (B0fac'*A*B0fac) ≈ 4π/3
	end



end


@testset "serial vs threaded" begin

    N = 10
    m = -N:N

    u = Inviscid(N; m)
    b = Insulating(N; m)

    B0s = [BasisElement(Basis{Insulating}, Poloidal, (1,0,1), 1.0), BasisElement(Basis{Insulating}, Toroidal, (1,0,1), 1.0)]
	U0s = [BasisElement(Basis{Inviscid}, Poloidal, (1,0,1), 1.0), BasisElement(Basis{Inviscid}, Toroidal, (1,0,1), 1.0)]

	@test Limace.diffusion(u) ≈ Limace.diffusion_threaded(u)
	@test Limace.diffusion(b) ≈ Limace.diffusion_threaded(b; external=true)

	@test Limace.inertial(u) ≈ Limace.inertial_threaded(u)
	@test Limace.inertial(b) ≈ Limace.inertial_threaded(b; external=true)

	@test Limace.coriolis(u) ≈ Limace.coriolis_threaded(u)


	for B0 in B0s
		RHSl = Limace.lorentz(u,b,B0)
		RHSi = Limace.induction(b,u,B0)
		RHSlt = Limace.lorentz_threaded(u,b,B0)
		RHSit = Limace.induction_threaded(b,u,B0)
		@test RHSl ≈ RHSlt
		@test RHSi ≈ RHSit
	end
	for U0 in U0s
		RHSl = Limace.lorentz(u,u,U0)
		RHSi = Limace.induction(b,U0,b)
		RHSlt = Limace.lorentz_threaded(u,u,U0)
		RHSit = Limace.induction_threaded(b,U0,b)
		@test RHSl ≈ RHSlt
		@test RHSi ≈ RHSit
	end

end