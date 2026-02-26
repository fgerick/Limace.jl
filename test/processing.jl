@testset "Processing: energies Malkus" begin

    N = 10
    Le = 1e-2

	for m in vcat(-N:-1, 1:N)
		u = Inviscid(N; m)
		b = PerfectlyConducting(N; m)

		B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15))

		bases = [u,b]


		forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u, b, B0), Limace.InductionB0(b,u,B0)]

		problem = LimaceProblem(bases, forcings)
		Limace.assemble!(problem)

		Limace.solve!(problem)

		ekin, emag = Limace.Processing.energies(problem)
		ekmratio = ekin./emag

		λ = problem.sol.values

		ekm_exact = @. abs(λ)^2/m^2 # Gerick (2020), eq. 3.16 without 8/15

		@test ekmratio ≈ ekm_exact
	end

end

@testset "Processing: misc" begin
	N0 = 10
	N1 = 5

    m = 0

	Le = 1e-3
	Lu = 2 / Le

    u = Inviscid(N0; m)
    b = Insulating(N0; m)

    B0 = BasisElement(b, Poloidal, (1,0,1))

    bases = [u,b]
    forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u,b,B0), Limace.InductionB0(b,u,B0), Limace.Diffusion(b, 1/Lu)]

    problem = LimaceProblem(bases, forcings)
    Limace.assemble!(problem; threads=true)

	a = problem.RHS

	u1 = Inviscid(N1; m)
	b1 = Insulating(N1; m)

    bases1 = [u1,b1]
    forcings1 = [Limace.Inertial(u1), Limace.Inertial(b1), Limace.Coriolis(u1, 1/Le), Limace.Lorentz(u1,b1,B0), Limace.InductionB0(b1,u1,B0), Limace.Diffusion(b1, 1/Lu)]

    problem1 = LimaceProblem(bases1, forcings1)
    Limace.assemble!(problem1; threads=true)

	a1 = problem1.RHS

	a_reduced = Limace.Processing.reducematrix(a, u, b, u1, b1)

	@test a1 ≈ a_reduced
end

@testset "Processing: spectra in l" begin

	N=4
	u = Inviscid(N; m=0)
	b = PerfectlyConducting(N; m=0)
	bases=[u,b]
	p = LimaceProblem(bases, [])
	nmat = size(p.LHS,1)
	λ = ones(ComplexF64,nmat)
	x = zeros(ComplexF64,nmat,nmat)

	nu = length(u)
	npu = Limace.Bases.np(u)
	npb = Limace.Bases.np(b)

	x[1,1] = 1.0 #u poloidal (1,0,0) = 1.0
	x[npu+2,1] = 2.0 #u toroidal (1,0,1) = 2.0
	x[nu+3,1] = 3.0 #b poloidal (2,0,0) = 3.0
	x[nu+npb+2,1] = 4.0 #b poloidal (1,0,1) = 4.0

	p.sol = GeneralizedEigen(λ,x)
	p.solved=true

	degs, specup, specut, specbp, specbt = Limace.Processing.spectrum(p; lmn=1)
	
	@test specup[2,1] ≈ 1/length(Limace.Bases._nrange_p(u,1))
	@test specut[2,1] ≈ 2^2/length(Limace.Bases._nrange_t(u,1))
	@test specbp[3,1] ≈ 3^2/length(Limace.Bases._nrange_p(b,2))
	@test specbt[2,1] ≈ 4^2/length(Limace.Bases._nrange_t(b,1))


end


@testset "Processing: eigenvector filter (functionality check only!)" begin

	N = 10
    m = 0

	Le = 1e-3
	Lu = 2 / Le

    u = Inviscid(N; m)
    b = Insulating(N; m)

    B0 = BasisElement(b, Poloidal, (1,0,1))

    bases = [u,b]
    forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u,b,B0), Limace.InductionB0(b,u,B0), Limace.Diffusion(b, 1/Lu)]

    problem = LimaceProblem(bases, forcings)
    Limace.assemble!(problem; threads=true)
	Limace.solve!(problem)

	nconverged = count(Limace.Processing.eigenvector_filter(problem; thresh=1e-2))
	@test nconverged < size(problem.LHS,1)
	nconverged = count(Limace.Processing.eigenvector_filter(problem; thresh=N))
	@test nconverged == size(problem.LHS,1)

	Ek = 1e-4
	u = Viscous(N; m=1)
	bases = [u]
    forcings = [Limace.Inertial(u), Limace.Coriolis(u), Limace.Diffusion(u, Ek)]

    problem = LimaceProblem(bases, forcings)
    Limace.assemble!(problem)
	Limace.solve!(problem)

	nconverged = count(Limace.Processing.eigenvector_filter(problem; thresh=1e-2))
	@test nconverged < size(problem.LHS,1)
	nconverged = count(Limace.Processing.eigenvector_filter(problem; thresh=N))
	@test nconverged == size(problem.LHS,1)
end