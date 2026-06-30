
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
	for b in (Limace.Inviscid, Limace.Insulating, Limace.Viscous, Limace.PerfectlyConducting)
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

    B0s = [BasisElement(b, Poloidal, (1,0,1), 1.0), BasisElement(b, Toroidal, (1,0,1), 1.0), BasisElement(b, Poloidal, (2,1,2), 1.0), BasisElement(b, Toroidal, (2,1,2), 1.0)]
	U0s = [BasisElement(u, Poloidal, (1,0,1), 1.0), BasisElement(u, Toroidal, (1,0,1), 1.0), BasisElement(u, Poloidal, (2,1,2), 1.0), BasisElement(u, Toroidal, (2,1,2), 1.0)]


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


@testset "Basis utils" begin
	N = 5
	u = Inviscid(N; m=0)
	nu = length(u)
	up100 = BasisElement(u, Poloidal, (1,0,0))
	ut100 = BasisElement(u, Toroidal, (1,0,0))

	@test size(u) == (nu,)
	@test u[1] == up100
	@test 2*up100 == BasisElement(u, Poloidal, (1,0,0), 2.0)
	@test -up100 == BasisElement(u, Poloidal, (1,0,0), -1.0)
	@test Limace.Bases.helmholtz(up100) == Poloidal
	@test Limace.Bases.helmholtz(ut100) == Toroidal 
	@test Limace.Bases.s(up100, u.V, 0.9) == Limace.Bases.s(u, 1,0,0, 0.9)
	@test Limace.Bases.t(ut100, u.V, 0.9) == Limace.Bases.t(u, 1,0,0, 0.9)


	# check missing implementations for custom basis
	struct TestB; end
	b = Limace.Basis{TestB, Limace.Sphere}(; N=1)
	@test_throws MethodError Limace.Bases.lpmax(b)
	@test_throws MethodError Limace.Bases.ltmax(b)
	@test_throws MethodError Limace.Bases.s(b, 1,0,0, 1.0)
	@test_throws MethodError Limace.Bases.t(b, 1,0,0, 1.0)
	@test_throws MethodError Limace.Bases.bcs_p(b)
	@test_throws MethodError Limace.Bases.bcs_t(b)
end

@testset "eigs" begin
    N = 10
    b = Limace.Inviscid(N)
    RHS = Limace.coriolis(b)
	max_eval = first(first(Limace.EigenSolve.eigs(RHS; nev=1)))
	max_eval2 = first(first(Limace.EigenSolve.eigs(RHS, sparse(1.0*I(size(RHS,1))); nev=1)))
	max_eval_dense = maximum(abs, eigvals(Matrix(RHS)))
	@test abs(max_eval) ≤ 2.0
	@test abs(max_eval) ≈ abs(max_eval_dense)
	@test abs(max_eval) ≈ abs(max_eval2)


end

@testset "poly derivatives" begin
	f = sin
	x = 0.1
	_f, _df = Limace.Poly.derivatives01(f,x)
	@test _f ≈ sin(x)
	@test _df ≈ cos(x)

	_f, _df, _d2f = Limace.Poly.derivatives012(f,x)
	@test _f ≈ sin(x)
	@test _df ≈ cos(x)
	@test _d2f ≈ -sin(x)

	_f, _df, _d2f, _d3f = Limace.Poly.derivatives0123(f,x)
	@test _f ≈ sin(x)
	@test _df ≈ cos(x)
	@test _d2f ≈ -sin(x)
	@test _d3f ≈ -cos(x)

end