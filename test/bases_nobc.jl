
@testset "Inviscid inertial modes BC explicit" begin

    function zhang(m, N) 
        sm = sign(m)
        m = abs(m)
        return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
    end

    N = 7
    u = Limace.Inviscid(N)
    RHS = Limace.coriolis(u)
    evals = eigvals(Matrix(RHS))

    u_nobc = Limace.InviscidNoBC(N+2)
    RHS_nobc = Limace.coriolis(u_nobc)
	Limace.boundarycondition!(RHS_nobc,u_nobc)
	LHS_nobc = Limace.inertial(u_nobc)
    evals_nobc = eigvals(Matrix(RHS_nobc),Matrix(LHS_nobc))

	for λ in evals
		@test any(x->norm(x-λ)≤ 1e-13,evals_nobc)
	end

end

@testset "Free decay modes BC explicit" begin
    #free decay modes damping
    λfd(l, k) = -besselj_zero(l - 1 / 2, k)[end]^2

    #get difference between λfd and targeted eigensolution
    function getlk1(Nmax, lmax, TB; external=true)
        b = TB(Nmax; m=0)
        LHS = Limace.inertial(b; external)
        RHS = Limace.diffusion(b; external)
		Limace.boundarycondition!(RHS,b)
        _evals = eigvals(Matrix(RHS),Matrix(LHS))
        for l = 1:lmax, k = 1:lmax÷2
            @test any(isapprox(λfd(l,k)), _evals)
        end
        return nothing
    end

    #testing:
    getlk1(20, 4, Limace.InsulatingNoBC; external=true)
    getlk1(20, 4, Limace.InsulatingNoBC; external=false)

    getlk1(20, 4, Limace.ThinWallBasis.ThinWall; external=true)
    getlk1(20, 4, Limace.ThinWallBasis.ThinWall; external=false)
end



@testset "Luo & Jackson 2022 mode NoBC" begin

    import Limace.Bases

	# struct LJ22; end
	# Limace.Bases.s(::Type{Basis{LJ22}}, V::Volume, l, m, n, r) = r^2 * (157 - 296r^2 + 143r^4) / (16 * sqrt(182 / 3))

	function compute(UT, BT; N=50)
		m = 0
		Le = 1e-4
		Lu = 2 / Le


		u = UT(N; m)
		b = BT(N; m)

		

		B0 = BasisElement(Basis{LJ22}, Poloidal, (2,0,1), 1.0)

		LHS = blockdiag(sparse(Limace.inertial(u),length(u),length(u)), sparse(Limace.inertial_threaded(b; external=false)))

		# @time begin
		RHSc = Limace.coriolis(u)/Le
		RHSl = Limace.lorentz_threaded(u,b,B0)
		RHSi = Limace.induction_threaded(b,u,B0; external=false)
		RHSd = Limace.diffusion_threaded(b;external=false)/Lu
		Limace.boundarycondition!(RHSd,b)
		# end

		RHS = [RHSc RHSl
			RHSi RHSd]

		# RHS = rhs(N, m; Ω = 2 / Le, η = 1 / Lu, lmnb0, B0poloidal = true, smfb0 = lj22)

		target = -0.0066 - 1.033im
		# target = -0.042+0.66im
		evals, _ = eigstarget(RHS, LHS, target; nev = 1)
		return first(evals)
	end


    lj22_n350 = -0.0065952461 - 1.0335959942im

	e = compute(Inviscid, Limace.InsulatingNoBC; N=60)
	@test e ≈ lj22_n350 atol=2e-4

	e = compute(Inviscid, Limace.ThinWall; N=60)
	@test e ≈ lj22_n350 atol=2e-4


end
