using Limace.Discretization: discretize, spectospat

@testset "discretize S₁⁰" begin
	b = Insulating(10)

	B0 = BasisElement(b, Poloidal, (1,0,1), sqrt(4pi/3)*Limace.InsulatingBasis.unitspherenorm(1,1))

	function b0_y10(r,θ,ϕ)	
		return (sqrt(7/46)*(3r^2-5)*cos(θ), sqrt(7/46)*(5-6r^2)*sin(θ), 0)
	end

	nr,nθ,nϕ = 50,50,50
	r = range(0.1,1.0,length=nr)
	θ = range(1e-9,π,length=nθ)
	ϕ = range(0,2π,length=nϕ)

	B1 = [discretize(B0, r, θ, ϕ) for r in r, θ in θ, ϕ in ϕ]
	B2 = [b0_y10(r, θ, ϕ) for r in r, θ in θ, ϕ in ϕ]


	for i=1:3
		@test getindex.(B1,i) ≈ getindex.(B2,i)
	end

end

@testset "spectospat S₁⁰" begin
	u = Inviscid(10)
	b = Insulating(10)

	x = zeros(ComplexF64, length(u)+length(b))

	x[length(u)+Limace.Bases.lmn2k_p_dict(b)[(1,0,1)]] = sqrt(4pi/3)*Limace.InsulatingBasis.unitspherenorm(1,1)

	function b0_y10(r,θ,ϕ)	
		return (sqrt(7/46)*(3r^2-5)*cos(θ), sqrt(7/46)*(5-6r^2)*sin(θ), 0)
	end


	nr,nθ,nϕ = 20,30,40
	ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(x, u, b, nr, nθ, nϕ)

	B2 = [b0_y10(r, π/2 -θ, ϕ) for r in r, θ in θ, ϕ in ϕ]
	br2, bθ2, bϕ2 = [getindex.(B2,i) for i=1:3]

	br3,bθ3,bϕ3, r,θ,ϕ = spectospat(x[length(u)+1:end], b, nr, nθ, nϕ)

	@test br ≈ br2 ≈ br3
	@test bθ ≈ bθ2 ≈ bθ3
	@test bϕ ≈ bϕ2 ≈ bϕ3

	# for one radius
	r=1.0
	ur,uθ,uϕ, br,bθ,bϕ, θ,ϕ = spectospat(x, u, b, r, nθ, nϕ)

	B2 = [b0_y10(r, π/2 -θ, ϕ) for r in r, θ in θ, ϕ in ϕ]
	br2, bθ2, bϕ2 = [getindex.(B2,i) for i=1:3]

	br3,bθ3,bϕ3, θ,ϕ = spectospat(x[length(u)+1:end], b, r, nθ, nϕ)

	@test br[1,:,:] ≈ br2 ≈ br3[1,:,:]
	@test bθ[1,:,:] ≈ bθ2 ≈ bθ3[1,:,:]
	@test bϕ[1,:,:] ≈ bϕ2 ≈ bϕ3[1,:,:]
end

@testset "spectospat (SHTns) vs. discretize (no SHTns)" begin
	N = 6
    Le = 1e-1

    u = Inviscid(N)
    b = PerfectlyConducting(N)

    bases = [u,b]

    B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15))

    forcings = [Limace.Inertial(u), Limace.Inertial(b), Limace.Coriolis(u, 1/Le), Limace.Lorentz(u, b, B0), Limace.InductionB0(b,u,B0)]

    problem = LimaceProblem(bases, forcings)
    Limace.assemble!(problem)
	Limace.solve!(problem)

	x = problem.sol.vectors[:,1]
	nr,nθ,nϕ = 20,30,40
	ur,uθ,uϕ, br,bθ,bϕ, r,θ,ϕ = spectospat(x, u, b, nr, nθ, nϕ)

	ur2,uθ2,uϕ2, br2, bθ2, bϕ2 = discretize(x, u, b, r, π/2 .- θ, ϕ)

	@test br ≈ br2
	@test bθ ≈ bθ2
	@test bϕ ≈ bϕ2
	@test ur ≈ ur2
	@test uθ ≈ uθ2
	@test uϕ ≈ uϕ2

end

