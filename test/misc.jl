
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

@testset "bases access functions" begin
	
	using Limace.Bases: lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax

	N = 10
	for b in (Limace.Inviscid, Limace.Insulating, Limace.Viscous)
		basis = b(N)
		@test basis.N == N
		@test basis.V == Limace.Bases.Sphere()
		@test basis.m == -N:N

		lmnp = lmn_p(basis)
		lmnt = lmn_t(basis)
		@test length(lmnp) == Limace.np(basis)
		@test length(lmnt) == Limace.nt(basis)
		@test length(lmnp)+length(lmnt) == length(basis)
		dp = lmn2k_p_dict(basis)
		dt = lmn2k_t_dict(basis)
		for (k,lmn) in enumerate(lmnp)
			@test dp[lmn] == k
		end
		for (k,lmn) in enumerate(lmnt)
			@test dt[lmn] == k
		end

	end

end

@testset "serial vs threaded" begin

    N = 10
    m = -N:N

    u = Inviscid(N; m)
    b = Insulating(N; m)

    B0s = [BasisElement(b, Poloidal, (1,0,1), 1.0), BasisElement(b, Toroidal, (1,0,1), 1.0)]
	U0s = [BasisElement(u, Poloidal, (1,0,1), 1.0), BasisElement(u, Toroidal, (1,0,1), 1.0)]


	@test Limace.diffusion(u) ≈ Limace._diffusion(Val(true),u)
	@test Limace.diffusion(b) ≈ Limace._diffusion(Val(true),b; external=true)

	@test Limace.inertial(u) ≈ Limace._inertial(Val(true),u)
	@test Limace.inertial(b) ≈ Limace._inertial(Val(true),b; external=true)

	@test Limace.coriolis(u) ≈ Limace._coriolis(Val(true),u)


	for B0 in B0s
		RHSl = Limace.lorentz(u,b,B0; threads=false)
		RHSi = Limace.induction(b,u,B0; threads=false)
		RHSlt = Limace.lorentz(u,b,B0; threads=true)
		RHSit = Limace.induction(b,u,B0; threads=true)
		@test RHSl ≈ RHSlt
		@test RHSi ≈ RHSit
	end
	for U0 in U0s
		RHSl = Limace.lorentz(u,u,U0; threads=false)
		RHSi = Limace.induction(b,U0,b; threads=false)
		RHSlt = Limace.lorentz(u,u,U0; threads=true)
		RHSit = Limace.induction(b,U0,b; threads=true)
		@test RHSl ≈ RHSlt
		@test RHSi ≈ RHSit
	end

end



@testset "eigs" begin
    N = 20
    b = Limace.Inviscid(N)
    RHS = Limace.coriolis(b)
	max_eval = first(first(Limace.Eigen.eigs(RHS; nev=1)))
	@test abs(max_eval) ≤ 2.0
end
