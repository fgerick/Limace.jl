using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod


Limace.DiscretePart.__wiginit(50)

@testset "Inviscid inertial modes" begin

    zhang(m,N)=-2/(m+2)*(√(1+m*(m+2)/(N*(2N+2m+1)))-1)*im

    N = 7
    RHS = Limace.InviscidBasis.rhs(N,1:N)
    evals = eigvals(Matrix(RHS))

    @test any(evals.≈1im)
    @test any(evals.≈zhang(1,1))
    @test any(evals.≈zhang(2,1))
    @test any(evals.≈zhang(3,1))
end


@testset "Free decay modes" begin
    #free decay modes damping
    λfd(l,k) = -besselj_zero(l-1/2,k)[end]^2

    function eigstarget(A,B,target; kwargs...)
        P = lu(A-target*B)
        LO = LinearMap{ComplexF64}((y,x)->ldiv!(y,P,B*x),size(A,2))
        pschur,history = partialschur(LO; kwargs...)
        evals, u = partialeigen(pschur)
        λ = 1 ./evals .+ target
        return λ,u
    end

    #get difference between λfd and targeted eigensolution
    function getlk1(Nmax,lmax)
        N,m,ns=Nmax,0,1:Nmax
        LHS = Limace.InsulatingMFBasis.lhs(N,m; ns)
        RHS = Limace.InsulatingMFBasis.rhs_diffusion(N,m; ns)
        found = [(eigstarget(RHS,LHS,λfd(l,k)+1e-12,nev=1)[1][1],λfd(l,k)) for l in 1:lmax for k in 1:lmax÷2]
        return found
    end

    #testing:
    f = getlk1(15,4)
    for (f_num,f_λfd) in f
        @test f_num ≈ f_λfd
    end
    
end

@testset "t₀¹s₀² kinematic dynamo" begin
    
    using Limace.DiscretePart: s_mf, t_mf, s_chen, t_chen, t_mf1, s_mf1
        
    # Li et al. (2010) eq. (24) fixed:
    @inline t_10(l,m,n,r) = 8.107929179422066*r*(1-r^2)
	@inline s_20(l,m,n,r) = 1.193271237996972*r^2*(1-r^2)^2

    
    N =24
    m=0
    LHS = Limace.InsulatingMFBasis.lhs(N,m)
    RHS_diff = Limace.InsulatingMFBasis.rhs_diffusion(N,m)

    ind_upol = Limace.DiscretePart.rhs_induction_upol(N,m,(2,0,1); su=s_20, smf = s_mf, tmf = t_mf,  condition=true)
    ind_utor = Limace.DiscretePart.rhs_induction_utor(N,m,(1,0,1); tu=t_10, smf = s_mf, tmf = t_mf, condition=true)

    Rm = 10
    RHS_ind = Rm * (ind_upol+ind_utor)
    RHS = RHS_diff + RHS_ind

    λdyn = eigvals(inv(Matrix(LHS))*Matrix(RHS))
    #get eigenvalue of largest real part
    λ = first(λdyn[sortperm(real.(λdyn),rev=true)]) 
    
    # reference value (Livermore & Jackson (2005), Table 1)
    λ_lj2005 = -8.01600246
    @test λ ≈ λ_lj2005
end

@testset "solid body rotation kinematic dynamo" begin
    
    
    using Limace.DiscretePart: s_mf, t_mf, s_in, t_in

    DP = Limace.DiscretePart
    N =12
    m=1
    
    LHS = Limace.InsulatingMFBasis.lhs(N,m)
    RHS_diff = Limace.InsulatingMFBasis.rhs_diffusion(N,m)
    DP.__wiginit(N)
    ind_utor = DP.rhs_induction_utor(N,m,(1,0,0); tu=t_in, smf = s_mf, tmf = t_mf, condition=false)

    Rm = 10
    RHS_ind = Rm * ind_utor
    RHS = RHS_diff + RHS_ind

    λdyn = eigvals(inv(Matrix(LHS))*Matrix(RHS))
    Limace.DiscretePart.wig_temp_free()


   #solid body rotation should add a constant imaginary part ( = frequency) to all eigenvalues
   
    @test all( imag.(λdyn).≈imag(λdyn[1]))

end


Limace.DiscretePart.wig_temp_free()