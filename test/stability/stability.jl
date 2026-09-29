
using Limace
using Limace: Inertial, Coriolis, Buoyancy, Diffusion, Lorentz, InductionB0, InductionU0, ScalarAdvectionT0, ScalarAdvectionU0, Advection

using BracketingNonlinearSolve



function update_Ra!(problem, Ra)
     i_f = findfirst(f->typeof(f)<:Buoyancy, problem.forcings)
    problem.forcings[i_f].factor = Ra
    Limace.assemble!(problem)
end

function update_and_solve_problem!(problem, α, update_f!)::ComplexF64
	update_f!(problem, α)
	Limace.solve!(problem; method=:sparse, target=0.0+1e-6im, which=:LR, nev=1) #, mindim=20, maxdim=300) 	
	λ = problem.sol.values
    if (length(λ)>0)             
		return first(λ)
	else
		return -1.0*Inf
	end
end

function solve_for_parameters(params; T0factor=2√π)
	(; Ek, Pr, N, m, Ra_interval) = params

	T = Limace.Temperature(N; m)
	u = Limace.Viscous(N; m)
	T0 = BasisElement(Basis{FP83, Limace.Sphere}, Toroidal, (0,0,3), T0factor)

	bases = [u,T]
	forcings = [Inertial(u), Coriolis(u, 1/Ek), Buoyancy(u,T), Diffusion(u),
				Inertial(T,Pr), ScalarAdvectionT0(T,u,T0), Diffusion(T)]

	problem = LimaceProblem(bases, forcings)
	Limace.assemble!(problem; threads=true)

	_ωc = [0.0]

	function f!(Ra,p)
		λ = update_and_solve_problem!(problem, Ra, update_Ra!)
		_ωc[1] = imag(λ)
		return real(λ)
	end
	interval_problem = IntervalNonlinearProblem(f!, Ra_interval)

	Ra_c = solve(interval_problem, ITP()).u	
	ω_c = _ωc[1]*Ek
	return Ra_c, ω_c
end

# T₀ = -r²/2, so that ∇T₀ = -𝐫
struct FP83; end
Limace.Bases.t(::Type{Basis{FP83, Limace.Sphere}}, V::Limace.Sphere, l, m, n, r) = -1/2*r^2


# Table 1 in https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.119.094501
# factor 2 in gravity?
@testset "Kaplan et al. (2017)" begin

	
	parameters = [(; Ek=1e-5, Pr=0.1, m=11, Ra_interval=(1e5,1e8), N=60),
				  (; Ek=3e-6, Pr=0.03, m=12, Ra_interval=(1e5,1e8), N=80),
				  (; Ek=1e-6, Pr=0.01, m=11, Ra_interval=(1e5,1e9), N=100)]
	
	sols = solve_for_parameters.(parameters; T0factor = 4√π)

	Ra_cs = getindex.(sols,1)
	ω_cs = getindex.(sols,2)

	Ra_cs_ref = [8.440e6, 2.336e7, 5.475e7]
	ω_cs_ref = [-0.04024, -0.04275, -0.03895]


	@test Ra_cs ≈ Ra_cs_ref rtol=5e-3
	@test ω_cs ≈ ω_cs_ref rtol=3e-3
	#tolerances can be lower for higher resolutions, but enough for CI.
	# Ek = 1e-7 case also not included due to need for high resolution.

end


@testset "Fearn & Proctor (1983), case a)" begin
	N = 20
	m = 2
	T = Limace.Temperature(N; m)
	u = Limace.Inviscid(N; m)
	b = Limace.Insulating(N; m)

	B0 = BasisElement(Basis{PerfectlyConducting, Limace.Sphere}, Toroidal, (1,0,0), 2sqrt(2pi/15))

	T0 = BasisElement(Basis{FP83, Limace.Sphere}, Toroidal, (0,0,3), 2√π)


	Λ = 0.5
	q = 1e-6

	bases = [u,b,T]
	forcings = [Coriolis(u, 1/2), Lorentz(u,b,B0,Λ), Buoyancy(u,T),
				Inertial(b,q), InductionB0(b,u,B0), Diffusion(b),
				Inertial(T), ScalarAdvectionT0(T,u,T0), Diffusion(T)]

	problem = LimaceProblem(bases, forcings)
	Limace.assemble!(problem)


	Ra_interval = (10.0,1e3)

	_f! = (Ra,p)->real(update_and_solve_problem!(problem, Ra, update_Ra!))


	interval_problem = IntervalNonlinearProblem(_f!, Ra_interval)

	Ra_c = solve(interval_problem, ITP()).u

	@test isapprox(115.4, Ra_c, atol=1e-1) #(3.6) in Fearn & Proctor (1983)

	update_Ra!(problem, Ra_c)
	Limace.solve!(problem; method=:sparse, target=0.0+0.0im, which=:LR, nev=1) 	
	λ = first(problem.sol.values)
	@test isapprox(43.14, abs(imag(λ)), atol = 1e-1) #(3.6) in Fearn & Proctor (1983)


end

# @testset "Fearn & Proctor (1983), case b)" begin
# 	N = 20
# 	m = 2
# 	T = Limace.Temperature(N; m)
# 	u = Limace.Inviscid(N; m)
# 	b = Limace.Insulating(N; m)

# 	B0 = [BasisElement(u, Toroidal, (2,0,0), 32/9*√(2π/105)),
# 		  BasisElement(u, Toroidal, (2,0,1), -32/9*√(2π/165))
# 		  ]	

# 	T0 = BasisElement(Basis{FP83, Limace.Sphere}, Toroidal, (0,0,3), 2√π)


# 	Λ = 1.0
# 	q = 1e-6

# 	bases = [u,b,T]
# 	forcings = [Coriolis(u, 1/2), Lorentz(u,b,B0,Λ), Buoyancy(u,T),
# 				Inertial(b,q), InductionB0(b,u,B0), Diffusion(b),
# 				Inertial(T), ScalarAdvectionT0(T,u,T0), Diffusion(T)]

# 	problem = LimaceProblem(bases, forcings)
# 	Limace.assemble!(problem)


# 	Ra_interval = (100.0,400.0)

# 	_f! = (Ra,p)->real(update_and_solve_problem!(problem, Ra))


# 	interval_problem = IntervalNonlinearProblem(_f!, Ra_interval)

# 	Ra_c = solve(interval_problem, ITP()).u

# 	@test isapprox(319.3, Ra_c, atol=1e-1) #Caption of Figure 3 in Fearn & Proctor (1983)

# 	update_Ra!(problem, Ra_c)
# 	Limace.solve!(problem; method=:sparse, target=0.0+0.0im, which=:LR, nev=1) 	
# 	λ = first(problem.sol.values)
# 	@test isapprox(36.54, abs(imag(λ)), atol = 1e-1) #Caption of Figure 3 in Fearn & Proctor (1983)


# end

@testset "Kore, onset viscous hydro, m=1,2,3" begin

	Ta = 1e9
	Ek = √(4/Ta)
	Pr = 1.0
	N = 50
	ms = 1:3
	Ra_interval = (1e6, 1e7)
	sols = [solve_for_parameters((; N, Ek, Pr, m, Ra_interval)) for m in 1:3]

	Ra_cs = getindex.(sols,1)
	ω_cs = getindex.(sols,2)


	#computed using kore (Ra_c, m, ω_c)
	kore_ref = [8.25010970e+06 1 -5.70072886e-03
				6.61663621e+06 2 -9.24957284e-03
				5.79607663e+06 3 -1.18798295e-02]

	Ra_cs_ref = kore_ref[:,1]
	ω_cs_ref = kore_ref[:,3]


	@test Ra_cs ≈ Ra_cs_ref rtol=1e-5
	@test ω_cs ≈ ω_cs_ref rtol=1e-5
	#tolerances can be lower for higher resolutions, but enough for CI.

end


#Table 1 in Maffei et al. (2024), https://doi.org/10.1093/gji/ggae294
@testset "Maffei et al. (2024)" begin
	parameters = [(; Ek=5e-4, Pr=0.1, m=3, Ra_interval=(1e4,1e7), N=30),
					(; Ek=1e-4, Pr=0.1, m=5, Ra_interval=(1e4,1e7), N=40),
					(; Ek=5e-5, Pr=0.1, m=6, Ra_interval=(1e5,1e7), N=60),
					(; Ek=1e-5, Pr=0.1, m=11, Ra_interval=(1e6,1e8), N=80),
					(; Ek=5e-6, Pr=0.1, m=13, Ra_interval=(1e6,1e8), N=100),
					(; Ek=1e-6, Pr=0.1, m=23, Ra_interval=(1e7,1e9), N=120),
					(; Ek=5e-4, Pr=1.0, m=4, Ra_interval=(1e5,1e7), N=30),
					(; Ek=5e-5, Pr=1.0, m=9, Ra_interval=(1e6,1e7), N=50)]

	sols = solve_for_parameters.(parameters)
	Ra_cs = getindex.(sols,1)
	ω_cs = getindex.(sols,2)

	Ra_cs_ref = [1.589e5, 9.642e5, 2.304e6, 1.688e7, 4.072e7, 3.245e8, 3.534e5, 6.308e6]
	ω_cs_ref = [-0.1213, -0.07855, -0.06354, -0.04025, -0.03145, -0.01909, -0.02878, -0.01812]


	for (Ra_c, Ra_c_ref) in zip(Ra_cs, Ra_cs_ref)
		@test Ra_c ≈ Ra_c_ref rtol=1e-3
	end
	for (ω_c, ω_c_ref) in zip(ω_cs, ω_cs_ref)
		@test ω_c ≈ ω_c_ref rtol=1e-3
	end
end

