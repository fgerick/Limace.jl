module MHDProblem

using Limace
using SparseArrays
using Wigxjpf

VB = Limace.InviscidBasis
BB = Limace.InsulatingMFBasis
DP = Limace.DiscretePart 

function lhs(N,m; ns = 0, Ω::T = 1.0, η::T = 1.0) where T
    LHSu = VB.lhs(N,m; ns, Ω)
    LHSb = BB.lhs(N,m; ns, Ω)

    LHS = blockdiag(LHSu,LHSb)
    return LHS
end

function rhs(N,m; ns = 0, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal = true, lmnb0::NTuple{3,Int} = (1,0,1), thresh = sqrt(eps()), kwargs...) where T
    

    @time "Coriolis" RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    # nu = size(RHSc,1)
    # @time "Viscous" if ν != 0
    #     RHSv = VB.rhs_viscosity(N,m; ns, ν)
    # else
    #     RHSv = spzeros(Complex{T}, nu, nu)
    # end

    RHSuu = RHSc #+ RHSv


    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

    if B0poloidal
        @time "Lorentz" RHSub = DP.rhs_lorentz_bpol(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_bpol(N,m, lmnb0; ns, η, thresh, kwargs...)
    else
        @time "Lorentz" RHSub = DP.rhs_lorentz_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

function lhs_cond(N,m; ns = 0, Ω::T = 1.0, η::T = 1.0) where T
    LHSu = VB.lhs(N,m; ns, Ω)
    LHSb = VB.lhs(N,m; ns, Ω)

    LHS = blockdiag(LHSu,LHSb)
    return LHS
end

function rhs_cond(N,m; ns = 0, Ω::T = 2.0, B0poloidal = false, B0fac=1.0, lmnb0::NTuple{3,Int} = (1,0,1), thresh = sqrt(eps()), kwargs...) where T
    

    RHSc = VB.rhs_coriolis(N,m; ns, Ω)

    RHSuu = RHSc #+ RHSv


    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

    if B0poloidal
        RHSub = DP.rhs_lorentz_bpol_cond(N,m, lmnb0; ns, η, thresh, kwargs...)*B0fac
        RHSbu = DP.rhs_induction_bpol_cond(N,m, lmnb0; ns, η, thresh, kwargs...)*B0fac
    else
        RHSub = DP.rhs_lorentz_btor_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
        RHSbu = DP.rhs_induction_btor_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

    RHS = [RHSuu RHSub
           RHSbu spzeros(size(RHSbu,1),size(RHSub,2))]
    
    return RHS
end
end