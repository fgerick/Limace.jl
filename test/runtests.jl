using Limace
using Test
using LinearAlgebra

zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im

@testset "Limace.jl" begin
    N = 7
    RHS = Limace.InviscidBasis.rhs(N,1:N)
    evals = eigvals(Matrix(RHS))

    @test any(evals.≈1im)
    @test any(evals.≈zhang(1,1))
    @test any(evals.≈zhang(2,1))
    @test any(evals.≈zhang(3,1))
end
