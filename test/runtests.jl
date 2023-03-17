using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod


Limace.DiscretePart.__wiginit(100)

@testset "Inviscid inertial modes" begin

    zhang(m, N) = -2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im

    N = 7
    RHS = Limace.InviscidBasis.rhs_coriolis(N, -N:N)
    evals = eigvals(Matrix(RHS))

    @test any(evals .≈ 1im)
    @test any(evals .≈ zhang(1, 1))
    @test any(evals .≈ zhang(2, 1))
    @test any(evals .≈ zhang(3, 1))
end


@testset "Free decay modes" begin
    #free decay modes damping
    λfd(l, k) = -besselj_zero(l - 1 / 2, k)[end]^2

    function eigstarget(A, B, target; kwargs...)
        P = lu(A - target * B)
        LO = LinearMap{ComplexF64}((y, x) -> ldiv!(y, P, B * x), size(A, 2))
        pschur, history = partialschur(LO; kwargs...)
        evals, u = partialeigen(pschur)
        λ = 1 ./ evals .+ target
        return λ, u
    end

    #get difference between λfd and targeted eigensolution
    function getlk1(Nmax, lmax)
        N, m, ns = Nmax, 0, 1:Nmax
        LHS = Limace.InsulatingMFBasis.lhs(N, m; ns)
        RHS = Limace.InsulatingMFBasis.rhs_diffusion(N, m; ns)
        found = [
            (eigstarget(RHS, LHS, λfd(l, k) + 1e-12, nev = 1)[1][1], λfd(l, k)) for
            l = 1:lmax for k = 1:lmax÷2
        ]
        return found
    end

    #testing:
    f = getlk1(15, 4)
    for (f_num, f_λfd) in f
        @test f_num ≈ f_λfd
    end

end



# @testset "t₀¹s₀¹ kinematic dynamo" begin

#     using Limace.DiscretePart: s_mf, t_mf, s_chen, t_chen, t_mf1, s_mf1

#     # Li et al. (2010) eq. (23):
#     @inline t_10(l,m,n,r) = r^2*sin(π*r)*√(4π/3)
# 	@inline s_10(l,m,n,r) = 0.17*t_10(l,m,n,r)


#     N =20
#     m=1
#     LHS = Limace.InsulatingMFBasis.lhs(N,m)
#     RHS_diff = Limace.InsulatingMFBasis.rhs_diffusion(N,m)

#     ind_upol = Limace.DiscretePart.rhs_induction_upol(N,m,(1,0,10); su=s_10, smf = s_mf, tmf = t_mf,  condition=false)
#     ind_utor = Limace.DiscretePart.rhs_induction_utor(N,m,(1,0,10); tu=t_10, smf = s_mf, tmf = t_mf, condition=false)

#     Rm = 160
#     RHS_ind = Rm * (ind_upol+ind_utor)
#     RHS = RHS_diff + RHS_ind

#     λdyn = eigvals(inv(Matrix(LHS))*Matrix(RHS))
#     #get eigenvalue of largest real part
#     λ = first(λdyn[sortperm(real.(λdyn),rev=true)]) 

#     # reference value (Li et al. (2010), Table 1)
#     λ_li2010 = 0.3131515894 - 34.8435927234im
#     @test λ ≈ λ_li2010
# end

@testset "solid body rotation kinematic dynamo" begin


    using Limace.DiscretePart: s_mf, t_mf, s_in, t_in

    DP = Limace.DiscretePart
    N = 12
    m = 2

    LHS = Limace.InsulatingMFBasis.lhs(N, m)
    RHS_diff = Limace.InsulatingMFBasis.rhs_diffusion(N, m)
    # DP.__wiginit(N)
    ind_utor = DP.rhs_induction_utor(
        N,
        m,
        (1, 0, 0);
        tu = t_in,
        smf = s_mf,
        tmf = t_mf,
        condition = false,
    )

    Rm = 10
    RHS_ind = Rm * ind_utor
    RHS = RHS_diff + RHS_ind

    λdyn = eigvals(inv(Matrix(LHS)) * Matrix(RHS))
    # Limace.DiscretePart.wig_temp_free()


    #solid body rotation should add a constant imaginary part ( = frequency) to all eigenvalues

    @test all(imag.(λdyn) .≈ imag(λdyn[1]))

end

@testset "t₀¹s₀² kinematic dynamo" begin

    using Limace.DiscretePart: s_mf, t_mf, s_chen, t_chen, t_mf1, s_mf1

    # Li et al. (2010) eq. (24) fixed:
    @inline t_10(l, m, n, r) = 8.107929179422066 * r * (1 - r^2)
    @inline s_20(l, m, n, r) = 1.193271237996972 * r^2 * (1 - r^2)^2

    function assemble(N)
        # N =45
        m = 0
        LHS = Limace.InsulatingMFBasis.lhs(N, m)
        RHS_diff = Limace.InsulatingMFBasis.rhs_diffusion(N, m)

        ind_upol = Limace.DiscretePart.rhs_induction_upol(
            N,
            m,
            (2, 0, 1);
            su = s_20,
            smf = s_mf,
            tmf = t_mf,
            condition = true,
        )
        ind_utor = Limace.DiscretePart.rhs_induction_utor(
            N,
            m,
            (1, 0, 1);
            tu = t_10,
            smf = s_mf,
            tmf = t_mf,
            condition = true,
        )
        return LHS, RHS_diff, ind_upol, ind_utor
    end

    function _RHS(Rm, RHS_diff, ind_upol, ind_utor)
        RHS_ind = Rm * (ind_upol + ind_utor)
        RHS = RHS_diff + RHS_ind
        return RHS
    end

    function getlargestEV(Rm, LHS, RHS_diff, ind_upol, ind_utor)
        RHS = _RHS(Rm, RHS_diff, ind_upol, ind_utor)

        λdyn = eigvals(inv(Matrix(LHS)) * Matrix(RHS))
        #get eigenvalue of largest real part
        λ = first(λdyn[sortperm(real.(λdyn), rev = true)])
        return λ
    end

    LHS, RHS_diff, ind_upol, ind_utor = assemble(45)
    # reference values (Livermore & Jackson (2005), Table 1,2)
    refvals = [(10, -8.01600246), (100, -6.92884871)]
    for (Rm, λ_lj2005) in refvals
        @test getlargestEV(Rm, LHS, RHS_diff, ind_upol, ind_utor) ≈ λ_lj2005
    end


end

@testset "Malkus modes" begin

    zhang(m, N) = -2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im

    # Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
    slow(m, N, Le, λ = imag(zhang(m, N))) =
        im * λ / 2Le * (1 - √(1 + 4Le^2 * m * (m - λ) / λ^2))
    fast(m, N, Le, λ = imag(zhang(m, N))) =
        im * λ / 2Le * (1 + √(1 + 4Le^2 * m * (m - λ) / λ^2))

    using Limace.MHDProblem: rhs_cond, lhs_cond
    N = 5
    m = 0:5
    Le = 1e-4

    lmnb0 = (1, 0, 0)

    LHS = lhs_cond(N, m)
    RHS = rhs_cond(
        N,
        m;
        Ω = 2 / Le,
        lmnb0,
        B0poloidal = false,
        B0fac = sqrt(4pi) / (sqrt(3) * sqrt(5 / 2)),
    )

    evals = eigvals(inv(Matrix(LHS)) * RHS)

    @test any(evals .≈ 1.0 * im / Le)
    for m = 1:(N-1)
        @test any(evals .≈ slow(m, 1, Le))
        @test any(evals .≈ fast(m, 1, Le))
    end

end

@testset "m=1 linear B₀ Mire.jl" begin
    #benchmark nonzonal m=1 background field 
    #against Mire.jl solution for B₀ = [0,1,0] × [x,y,z]

    #Mire.jl code:
    #
    #prob = MHDProblem(3, Sphere(), [0,0,1.0], 1e-3, [0,1,0] × Mire.𝐫, IversBasis, IversBasis); 
    #assemble!(prob);
    #λ_mire,u_mire = eigen(inv(Matrix(prob.LHS))*prob.RHS);

    λ_mire = ComplexF64[
        -0.005095583353799794+0.0im,
        -0.00447202108007672+0.0im,
        -2.9796487116893485e-9+0.0im,
        -8.080392367089447e-11+0.0im,
        -5.950795411990839e-14-1708.0249665649847im,
        -5.950795411990839e-14+1708.0249665649847im,
        -5.906386491005833e-14-611.9886647220552im,
        -5.906386491005833e-14+611.9886647220552im,
        -5.684341886080802e-14-1509.9412787288688im,
        -5.684341886080802e-14+1509.9412787288688im,
        -2.3276519600656798e-14-0.0018617314109947115im,
        -2.3276519600656798e-14+0.0018617314109947115im,
        -1.0658141036401503e-14-176.61194539178334im,
        -1.0658141036401503e-14+176.61194539178334im,
        -9.020562075079397e-17-1309.309123548131im,
        -9.020562075079397e-17+1309.309123548131im,
        -1.3877787807814457e-17-231.931995092941im,
        -1.3877787807814457e-17+231.931995092941im,
        0.0-820.0114646508597im,
        0.0+0.0im,
        0.0+0.0im,
        0.0+0.0im,
        0.0+0.0im,
        0.0+0.0im,
        0.0+0.0im,
        0.0+820.0114646508597im,
        1.3877787807814457e-17-1231.9269950324663im,
        1.3877787807814457e-17+1231.9269950324663im,
        4.99437628008817e-17+0.0im,
        2.0304667235765406e-15-0.0019302965014471438im,
        2.0304667235765406e-15+0.0019302965014471438im,
        2.227493363371469e-15-0.0065528310972463415im,
        2.227493363371469e-15+0.0065528310972463415im,
        2.3105161447295464e-15-1.5275119964485748im,
        2.3105161447295464e-15+1.5275119964485748im,
        2.635017854954458e-15-0.00199999337503616im,
        2.635017854954458e-15+0.00199999337503616im,
        6.727717051788765e-15-1.1585515372080849e-8im,
        6.727717051788765e-15+1.1585515372080849e-8im,
        1.8437122253356578e-14+0.0im,
        3.4638958368304884e-14-1000.0000000000002im,
        3.4638958368304884e-14+1000.0000000000002im,
        3.552713678800501e-14-500.0024999163528im,
        3.552713678800501e-14+500.0024999163528im,
        8.526512829121202e-14-666.6676666560261im,
        8.526512829121202e-14+666.6676666560261im,
        1.4210854715202004e-13-894.4283090383414im,
        1.4210854715202004e-13+894.4283090383414im,
        8.080512003078941e-11+0.0im,
        2.979634728476476e-9+0.0im,
        0.004472021080068012+0.0im,
        0.005095583353808171+0.0im,
    ]

    N = 3
    m = -N:N
    Le = 1e-3


    LHS = lhs_cond(N, m)
    RHS = rhs_cond(
        N,
        m;
        Ω = 2 / Le,
        lmnb0 = (1, 1, 0),
        B0poloidal = false,
        B0fac = sqrt(16pi / 15),
    )

    λ, u = eigen(Matrix(RHS))

    @test sort(imag.(λ)) ≈ sort(imag.(λ_mire))
end
@testset "Luo & Jackson 2022 mode" begin

    function eigstarget(A, B, target; kwargs...)
        P = lu(A - target * B)
        LO = LinearMap{ComplexF64}((y, x) -> ldiv!(y, P, B * x), size(A, 2))
        pschur, history = partialschur(LO; kwargs...)
        evals, u = partialeigen(pschur)
        λ = 1 ./ evals .+ target
        return λ, u
    end

    using Limace.MHDProblem: rhs, lhs
    N = 50
    m = 0
    Le = 1e-4
    Lu = 2 / Le

    lmnb0 = (2, 0, 2)
    lj22(l, m, n, r) = r^2 * (157 - 296r^2 + 143r^4) / (16 * sqrt(182 / 3))
    # lj22(l,m,n,r) = r^2*(157-296r^2+143r^4)*5/14*sqrt(3/182)

    # lmnb0 = (1,0,1)
    # lj22(l,m,n,r) = r*(5-3r^2)*sqrt(7/46)/2

    LHS = lhs(N, m)
    RHS = rhs(N, m; Ω = 2 / Le, η = 1 / Lu, lmnb0, B0poloidal = true, smfb0 = lj22)

    target = -0.0066 - 1.033im
    # target = -0.042+0.66im
    evals, evecs = eigstarget(RHS, LHS, target; nev = 4)

    lj22_n350 = -0.0065952461 - 1.0335959942im

    @test any(isapprox.(evals, lj22_n350, atol = 1e-7))


end

Limace.DiscretePart.wig_temp_free()