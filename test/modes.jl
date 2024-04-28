
function eigstarget(A, B, target; kwargs...)
    P = lu(A - target * B)
    LO = LinearMap{ComplexF64}((y, x) -> ldiv!(y, P, B * x), size(A, 2))
    pschur, history = partialschur(LO; kwargs...)
    evals, u = partialeigen(pschur)
    λ = 1 ./ evals .+ target
    return λ, u
end


@testset "Inviscid inertial modes" begin

    # zhang(m, N) = -2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
    function zhang(m, N) 
        sm = sign(m)
        m = abs(m)
        return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
    end

    N = 7
    b = Limace.Inviscid(N)
    RHS = Limace.coriolis(b)
    evals = eigvals(Matrix(RHS))

    @test any(evals .≈ 1im)
    @test any(evals .≈ zhang(1, 1))
    @test any(evals .≈ zhang(2, 1))
    @test any(evals .≈ zhang(3, 1))

    for m in -(N-1):(N-1)
        b = Limace.Inviscid(N; m)
        RHS = Limace.coriolis(b)
        evals = eigvals(Matrix(RHS))
        m == 1 && @test any(evals .≈ 1im)
        m == -1 && @test any(evals .≈ -1im)
        @test any(evals .≈ zhang(m, 1))
    end
end

@testset "Inviscid inertial modes, rotating frame vs. u0 = Ωs𝐞ᵩ" begin
    function assembleU0uniHydro(N,m; ns=false, Ω=2.0, Ω2 = 0.0)

        u = Inviscid(N; m)
        RHS = Limace.coriolis(u; Ω)
        if Ω2 != 0.0
            Limace.Poly.__wiginit(N)
            U0 = BasisElement(Basis{Inviscid}, Toroidal, (1,0,0), 2sqrt(2pi/15))
            RHSadv = -Ω2*Limace.lorentz(u, u, U0) #advection term is the same as Lorentz term
            RHS += RHSadv
        end
        return RHS
    end
    N = 3
    for m = -3:3
        λ = eigvals(Matrix(assembleU0uniHydro(N,m)))
        λu0 = eigvals(Matrix(assembleU0uniHydro(N,m; Ω=0.0, Ω2 = 1.0))) .+ m*im

        @test λ[sortperm(imag.(λ))] ≈ λu0[sortperm(imag.(λu0))] 
    end

end

@testset "Free decay modes" begin
    #free decay modes damping
    λfd(l, k) = -besselj_zero(l - 1 / 2, k)[end]^2

    #get difference between λfd and targeted eigensolution
    function getlk1(Nmax, lmax)
        b = Limace.Insulating(Nmax; n=1:Nmax)
        LHS = Limace.inertial(b)
        RHS = Limace.diffusion(b)
        found = [
            (first(first(eigstarget(RHS, LHS, λfd(l, k) + 1e-9, nev = 1))), λfd(l, k)) for
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

#     using Limace.DiscretePart: s_mf, t_mf, s_chen, t_chen

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


    # DP = Limace.DiscretePart
    N = 12
    m = 2
    b = Insulating(N; m)

    Limace.Poly.__wiginit(2N)

    LHS = Limace.inertial(b)
    RHS_diff = Limace.diffusion(b)
    U0 = BasisElement(Basis{Inviscid}, Toroidal, (1,0,0), 1.0)
    ind_utor = Limace.induction(b,U0,b)

    Rm = 10.0
    RHS_ind = Rm * ind_utor
    RHS = RHS_diff + RHS_ind

    λdyn = eigvals(LHS\Matrix(RHS))
    Limace.Poly.wig_temp_free()

    #solid body rotation should add a constant imaginary part ( = frequency) to all eigenvalues

    @test all(imag.(λdyn) .≈ imag(λdyn[1]))

end

@testset "t₀¹s₀² kinematic dynamo" begin

    import Limace.Bases


    struct Li2010; end
    # Li et al. (2010) eq. (24) fixed:
    @inline Limace.Bases.t(::Type{Basis{Li2010}}, l, m, n, r) = 8.107929179422066 * r * (1 - r^2) #t10
    @inline Limace.Bases.s(::Type{Basis{Li2010}},l, m, n, r) = 1.193271237996972 * r^2 * (1 - r^2)^2 #s20

    function assemble(N)
        # N =45
        Limace.Poly.__wiginit(2N)

        m = 0
        b = Insulating(N; m)
        LHS = Limace.inertial(b)
        RHS_diff = Limace.diffusion(b)

        U0_t10 = BasisElement(Basis{Li2010}, Toroidal, (1,0,1), 1.0)
        U0_s20 = BasisElement(Basis{Li2010}, Poloidal, (2,0,1), 1.0)

        RHS_induction_t10 = Limace.induction(b,U0_t10,b)
        RHS_induction_s20 = Limace.induction(b,U0_s20,b)
        Limace.Poly.wig_temp_free()
        return LHS, RHS_diff, RHS_induction_s20, RHS_induction_t10
    end

    function _RHS(Rm, RHS_diff, ind_upol, ind_utor)
        RHS_ind = Rm * (ind_upol + ind_utor)
        RHS = RHS_diff + RHS_ind
        return RHS
    end

    function getlargestEV(Rm, LHS, RHS_diff, ind_upol, ind_utor)
        RHS = _RHS(Rm, RHS_diff, ind_upol, ind_utor)

        λdyn = eigvals(LHS\Matrix(RHS))
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

    function zhang(m, N) 
        sm = sign(m)
        m = abs(m)
        return -sm*2 / (m + 2) * (√(1 + m * (m + 2) / (N * (2N + 2m + 1))) - 1) * im
    end

    # Malkus J. Fluid Mech. (1967), vol. 28, pp. 793-802, eq. (2.28)
    slow(m, N, Le, λ = imag(zhang(m, N))) =
        im * λ / 2Le * (1 - √(1 + 4Le^2 * m * (m - λ) / λ^2))
    fast(m, N, Le, λ = imag(zhang(m, N))) =
        im * λ / 2Le * (1 + √(1 + 4Le^2 * m * (m - λ) / λ^2))

    # using Limace.MHDProblem: rhs_cond_pre, lhs_cond
    N = 6
    Limace.Poly.__wiginit(2N)
    # m = -5:5
    Le = 1e-1

    u = Inviscid(N)
    b = PerfectlyConducting(N)


    B0 = BasisElement(b, Toroidal, (1,0,0), 2sqrt(2pi/15))

   RHS = [Limace.coriolis(u)/Le Limace.lorentz(u, b, B0)
          Limace.induction(b,u,B0) spzeros(length(b),length(u))];

    evals = eigvals(Matrix(RHS))

    @test any(evals .≈ 1.0 * im / Le)
    for m = vcat(-(N-1):-1, 1:(N-1))
        @test any(evals .≈ slow(m, 1, Le))
        @test any(evals .≈ fast(m, 1, Le))
    end

    for m = vcat(-(N-1):-1, 1:(N-1))
        u = Inviscid(N; m)
        b = PerfectlyConducting(N; m)
        RHS = [Limace.coriolis(u)/Le Limace.lorentz(u,b, B0)
        Limace.induction(b,u,B0) spzeros(length(b),length(u))];
        evals = eigvals(Matrix(RHS))    
        @test any(evals .≈ slow(m, 1, Le))
        @test any(evals .≈ fast(m, 1, Le))
    end

    Limace.Poly.wig_temp_free()
end

# @testset "Malkus modes, rotating frame vs. u0 = Ωs𝐞ᵩ" begin

#     function assembleU0uniMalkus(N,m; ns=false, Ω=2.0, Ω2 = 0.0)
#         r,wr = DP.rquad(N+5+5)
    
#         js_a1 = DP.jacobis_l(N,r,1.0)
#         js_a0 = DP.jacobis_l(N,r,0.0)

#         RHSuu = Limace.InviscidBasis.rhs_coriolis(N,m; ns, Ω)

#         DP.wig_table_init(2N, 9)
#         DP.wig_temp_init(2N)
#         B0fac = 2sqrt(2pi/15)
#         RHSub = DP.rhs_lorentz_btor_cond_pre(N,m, (1,0,0), r, wr, js_a1, js_a0; ns)*B0fac
#         RHSbu = DP.rhs_induction_btor_cond_pre(N,m, (1,0,0), r, wr, js_a1, js_a0; ns)*B0fac
#         RHSbb = spzeros(size(RHSbu,1),size(RHSub,2))
#         if Ω2 != 0.0
#             RHSuu += -Ω2*2sqrt(2pi/15)*DP.rhs_advection_utor_pre(N,m, (1,0,0), r,wr, js_a1,js_a0, conditions=true) 
#             RHSbb += Ω2*2sqrt(2pi/15)*DP.rhs_induction_utor_cond_pre(N,m, (1,0,0), r,wr, js_a1,js_a0; ns) 
#         end


#         DP.wig_temp_free()
#         RHS = [RHSuu RHSub
#                 RHSbu RHSbb]
#         return RHS
#     end
#     N = 3
#     Le = 1e-2
#     for m = -3:3
#         λ = eigvals(Matrix(assembleU0uniMalkus(N,m; Ω=2/Le)))
#         λu0 = eigvals(Matrix(assembleU0uniMalkus(N,m; Ω=0.0, Ω2 = 1/Le))) .+ m*im/Le
#         @test λ[sortperm(imag.(λ))] ≈ λu0[sortperm(imag.(λu0))]
#     end

# end

# @testset "m=1 linear B₀ Mire.jl" begin
#     #benchmark nonzonal m=1 background field 
#     #against Mire.jl solution for B₀ = [0,1,0] × [x,y,z]

#     #Mire.jl code:
#     #
#     #prob = MHDProblem(3, Sphere(), [0,0,1.0], 1e-3, [0,1,0] × Mire.𝐫, IversBasis, IversBasis); 
#     #assemble!(prob);
#     #λ_mire,u_mire = eigen(inv(Matrix(prob.LHS))*prob.RHS);

#     λ_mire = ComplexF64[
#         -0.005095583353799794+0.0im,
#         -0.00447202108007672+0.0im,
#         -2.9796487116893485e-9+0.0im,
#         -8.080392367089447e-11+0.0im,
#         -5.950795411990839e-14-1708.0249665649847im,
#         -5.950795411990839e-14+1708.0249665649847im,
#         -5.906386491005833e-14-611.9886647220552im,
#         -5.906386491005833e-14+611.9886647220552im,
#         -5.684341886080802e-14-1509.9412787288688im,
#         -5.684341886080802e-14+1509.9412787288688im,
#         -2.3276519600656798e-14-0.0018617314109947115im,
#         -2.3276519600656798e-14+0.0018617314109947115im,
#         -1.0658141036401503e-14-176.61194539178334im,
#         -1.0658141036401503e-14+176.61194539178334im,
#         -9.020562075079397e-17-1309.309123548131im,
#         -9.020562075079397e-17+1309.309123548131im,
#         -1.3877787807814457e-17-231.931995092941im,
#         -1.3877787807814457e-17+231.931995092941im,
#         0.0-820.0114646508597im,
#         0.0+0.0im,
#         0.0+0.0im,
#         0.0+0.0im,
#         0.0+0.0im,
#         0.0+0.0im,
#         0.0+0.0im,
#         0.0+820.0114646508597im,
#         1.3877787807814457e-17-1231.9269950324663im,
#         1.3877787807814457e-17+1231.9269950324663im,
#         4.99437628008817e-17+0.0im,
#         2.0304667235765406e-15-0.0019302965014471438im,
#         2.0304667235765406e-15+0.0019302965014471438im,
#         2.227493363371469e-15-0.0065528310972463415im,
#         2.227493363371469e-15+0.0065528310972463415im,
#         2.3105161447295464e-15-1.5275119964485748im,
#         2.3105161447295464e-15+1.5275119964485748im,
#         2.635017854954458e-15-0.00199999337503616im,
#         2.635017854954458e-15+0.00199999337503616im,
#         6.727717051788765e-15-1.1585515372080849e-8im,
#         6.727717051788765e-15+1.1585515372080849e-8im,
#         1.8437122253356578e-14+0.0im,
#         3.4638958368304884e-14-1000.0000000000002im,
#         3.4638958368304884e-14+1000.0000000000002im,
#         3.552713678800501e-14-500.0024999163528im,
#         3.552713678800501e-14+500.0024999163528im,
#         8.526512829121202e-14-666.6676666560261im,
#         8.526512829121202e-14+666.6676666560261im,
#         1.4210854715202004e-13-894.4283090383414im,
#         1.4210854715202004e-13+894.4283090383414im,
#         8.080512003078941e-11+0.0im,
#         2.979634728476476e-9+0.0im,
#         0.004472021080068012+0.0im,
#         0.005095583353808171+0.0im,
#     ]

#     N = 3
#     m = -N:N
#     Le = 1e-3


#     LHS = lhs_cond(N, m)
#     RHS = rhs_cond_pre(
#         N,
#         m;
#         Ω = 2 / Le,
#         lmnb0 = (1, 1, 0),
#         B0poloidal = false,
#         B0fac = sqrt(16pi / 15),
#     )

#     λ, u = eigen(Matrix(RHS))

#     @test sort(imag.(λ)) ≈ sort(imag.(λ_mire))
# end

@testset "Luo & Jackson 2022 mode" begin

    import Limace.Bases

    N = 50
    m = 0
    Le = 1e-4
    Lu = 2 / Le

    Limace.Poly.__wiginit(N)

    u = Inviscid(N; m)
    b = Insulating(N; m)

    struct LJ22; end

    Limace.Bases.s(::Type{Basis{LJ22}}, l, m, n, r) = r^2 * (157 - 296r^2 + 143r^4) / (16 * sqrt(182 / 3))

    B0 = BasisElement(Basis{LJ22}, Poloidal, (2,0,1), 1.0)

    LHS = SymTridiagonal(blockdiag(sparse(Limace.inertial(u),length(u),length(u)), sparse(Limace.inertial(b))))

    RHSc = Limace.coriolis(u)/Le
    RHSl = Limace.lorentz(u,b,B0)
    RHSi = Limace.induction(b,u,B0)
    RHSd = Limace.diffusion(b)/Lu

    RHS = [RHSc RHSl
           RHSi RHSd]
    # RHS = rhs(N, m; Ω = 2 / Le, η = 1 / Lu, lmnb0, B0poloidal = true, smfb0 = lj22)

    target = -0.0066 - 1.033im
    # target = -0.042+0.66im
    evals, evecs = eigstarget(RHS, LHS, target; nev = 1)

    lj22_n350 = -0.0065952461 - 1.0335959942im

    abs(lj22_n350-first(evals))
    @test any(isapprox.(evals, lj22_n350, atol = 1e-7))

    Limace.Poly.wig_temp_free()

end

@testset "Luo & Jackson 2022 dipole mode" begin

    N = 130
    m = 0


	Le = 1e-3
	Lu = 2 / Le

    u = Inviscid(N; m)
    b = Insulating(N; m)

    B0 = BasisElement(Basis{Insulating}, Poloidal, (1,0,1), sqrt(30/23))

    LHS = SymTridiagonal(blockdiag(sparse(Limace.inertial(u),length(u),length(u)), sparse(Limace.inertial(b))))

    Limace.Poly.__wiginit(N)
    RHSc = Limace.coriolis(u)/Le
    RHSl = Limace.lorentz(u,b,B0)
    RHSi = Limace.induction(b,u,B0)
    RHSd = Limace.diffusion(b)/Lu
    Limace.Poly.wig_temp_free()

    RHS = [RHSc RHSl
           RHSi RHSd]

	target = -0.041950864156977755 - 0.6599458208985812im
	evals, evecs = eigstarget(RHS, LHS, target; nev = 1)
	@test isapprox(first(evals), target, atol = 1e-6)

end


@testset "Luo, Marti & Jackson 2022 s₁⁰" begin


    N = 40
    m = 1

    Eη = 1e-9
    Le = 2√(Eη)
    Lu = 2 / Le

    lmnb0 = (1, 0, 1)
    B0fac = -4sqrt(pi/35)
  
    u = Inviscid(N; m)
    b = Insulating(N; m)
    B0 = BasisElement(Basis{Insulating}, Poloidal, lmnb0, B0fac)

    LHSu = sparse(Limace.inertial(u),length(u),length(u))*Eη
    LHSb = sparse(Limace.inertial(b))
    LHS = SymTridiagonal(blockdiag(LHSu,LHSb))

    Limace.Poly.__wiginit(N)
    RHSc = Limace.coriolis(u; Ω = 1.0)
    RHSl = Limace.lorentz(u,b,B0)
    RHSi = Limace.induction(b,u,B0)
    RHSd = Limace.diffusion(b)
    Limace.Poly.wig_temp_free()

    RHS = [RHSc RHSl
           RHSi RHSd]

    target = -287.9448432-115.2081087im

    evals, evecs = eigstarget(RHS, LHS, target; nev = 1)

    @test isapprox(first(evals), target, atol = 1e-4) #at N=40 we match the 1e-4 converged digits of N=120 of LMJ2022


end

@testset "Luo, Marti & Jackson 2022 t₁⁰" begin


    N = 70
    m = 3
    Eη = 1e-9
    Le = 2√(Eη)
    Lu = 2 / Le

    lmnb0 = (1, 0, 1)


    B0fac = -4sqrt(pi/35)
    
    u = Inviscid(N; m)
    b = Insulating(N; m)
    B0 = BasisElement(Basis{Insulating}, Toroidal, lmnb0, B0fac)

    LHSu = sparse(Limace.inertial(u),length(u),length(u))*Eη
    LHSb = sparse(Limace.inertial(b))
    LHS = SymTridiagonal(blockdiag(LHSu,LHSb))

    Limace.Poly.__wiginit(N)
    RHSc = Limace.coriolis(u; Ω = 1.0)
    RHSl = Limace.lorentz(u,b,B0)
    RHSi = Limace.induction(b,u,B0)
    RHSd = Limace.diffusion(b)
    Limace.Poly.wig_temp_free()

    RHS = [RHSc RHSl
           RHSi RHSd]


    target = -742.7652176+684.132152im
    evals, evecs = eigstarget(RHS, LHS, target; nev = 1)

    @test isapprox(first(evals), target, atol = 1e-4) #at N=70 we match the 1e-4 converged digits of N=120 of LMJ2022


end

@testset "Luo, Marti & Jackson 2022 t₁⁰s₁⁰" begin


    N = 70
    m = 5
    Eη = 1e-9
    Le = 2√(Eη)
    Lu = 2 / Le

    lmnb0 = (1, 0, 1)


    B0fact = -1/√2
    B0facp = -√(15/23)
   
    u = Inviscid(N; m)
    b = Insulating(N; m)
    B0t = BasisElement(Basis{Insulating}, Toroidal, lmnb0, B0fact)
    B0p = BasisElement(Basis{Insulating}, Poloidal, lmnb0, B0facp)

    LHSu = sparse(Limace.inertial(u),length(u),length(u))*Eη
    LHSb = sparse(Limace.inertial(b))
    LHS = SymTridiagonal(blockdiag(LHSu,LHSb))

    Limace.Poly.__wiginit(N)
    RHSc = Limace.coriolis(u; Ω = 1.0)
    RHSl = Limace.lorentz(u,b,B0t) + Limace.lorentz(u,b, B0p)
    RHSi = Limace.induction(b,u,B0t) + Limace.induction(b,u, B0p)
    RHSd = Limace.diffusion(b)
    Limace.Poly.wig_temp_free()

    RHS = [RHSc RHSl
    RHSi RHSd]

    target = -2134.0 + 2555.2im
    evals, evecs = eigstarget(RHS, LHS, target; nev = 1)

    @test isapprox(first(evals), target, atol = 1e-1) # only one significant digit available


end

# @testset "Distributed vs serial" begin
#     addprocs(4; exeflags=`--project=$(Base.active_project())`)
#     @everywhere begin
#         using Limace
#         using Limace.MHDProblem: rhs, rhs_dist
#         N = 10
#         m = -N:N
#     end
#     @test rhs_dist(N,m) ≈ rhs(N,m)
#     rmprocs(workers())
# end
