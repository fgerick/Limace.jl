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