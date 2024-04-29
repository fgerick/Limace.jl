
@testset "Inviscid inertial modes" begin

    function zhang(m, N) 
        sm = sign(m)
        m = abs(m)
        return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
    end

    N = 7
    u = Limace.Inviscid(N)
    RHS = Limace.coriolis(u)
    evals = eigvals(Matrix(RHS))

    u_nobc = Limace.InviscidNoBC(N+3)
    RHS_nobc = Limace.coriolis(u_nobc)
	LHS_nobc = Limace.inertial(u_nobc)
    evals_nobc = eigvals(Matrix(RHS_nobc),Matrix(LHS_nobc))

	for λ in evals
		@test any(x->norm(x-λ)≤ 1e-13,evals_nobc)
	end

end