
@testset "unit sphere norm of n=1 poloidal insulating MF fields" begin
	N = 20
	r,wr = DP.rquad(N)

	for l in 1:5, m in -l:l, n in 1:3
		@test DP._inertial_SS((l,m,n), (l,m,n), r,wr, DP.s_mf,DP.s_mf)*Limace.InsulatingMFBasis.unitspherenorm(l,n)^2 ≈ 1.0
	end
end

@testset "normalizing non-orthogonal background fields so that total (4π/3)⁻¹ ∫₀¹  B⋅B dV = 1" begin

	using Limace.InsulatingMFBasis: unitspherenorm, norm_B0fac!, inner_b0norm

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
		@test 4π/3*(B0fac'*A*B0fac) ≈ 1
	end



end