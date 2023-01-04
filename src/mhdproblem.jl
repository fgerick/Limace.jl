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
    nu = size(RHSc,1)
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
        @time "Induction" RHSbu = DP.rhs_induction_bpol(N,m, lmnb0; ns, η, thresh)
    else
        @time "Lorentz" RHSub = DP.rhs_lorentz_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_btor(N,m, lmnb0; ns, η, thresh)
    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

end