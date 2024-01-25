
@testset "unit sphere norm of n=1 poloidal insulating MF fields" begin
	N = 20
	r,wr = DP.rquad(N)

	for l in 1:5, m in -l:l, n in 1:3
		@test DP._inertial_SS((l,m,n), (l,m,n), r,wr, DP.s_mf,DP.s_mf)*Limace.InsulatingMFBasis.unitspherenorm(l,n)^2 ≈ 1.0
	end
end