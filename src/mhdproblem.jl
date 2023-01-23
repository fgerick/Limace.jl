module MHDProblem

using Limace
using SparseArrays
using Wigxjpf
using Distributed

VB = Limace.InviscidBasis
BB = Limace.InsulatingMFBasis
DP = Limace.DiscretePart 

function lhs(N,m; ns = 0, Ω::T = 1.0, η::T = 1.0) where T
    LHSu = VB.lhs(N,m; ns, Ω)
    LHSb = BB.lhs(N,m; ns, Ω)

    LHS = blockdiag(LHSu,LHSb)
    return LHS
end

function rhs(N,m; ns = 0, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal = true, lmnb0::NTuple{3,Int} = (1,0,1), thresh = 1000eps(), kwargs...) where T
    
    lb0,mb0,nb0 = lmnb0
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
        if mb0 != 0
            @time "Lorentz" RHSub += (-1)^mb0*DP.rhs_lorentz_bpol(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            @time "Induction" RHSbu += (-1)^mb0*DP.rhs_induction_bpol(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        @time "Lorentz" RHSub = DP.rhs_lorentz_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
        if mb0 != 0
            @time "Lorentz" RHSub += (-1)^mb0*DP.rhs_lorentz_btor(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            @time "Induction" RHSbu += (-1)^mb0*DP.rhs_induction_btor(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

function rhs_dist(N,m; ns = 0, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal = true, lmnb0::NTuple{3,Int} = (1,0,1), thresh = 1000eps(), kwargs...) where T
    
    lb0,mb0,nb0 = lmnb0
    @time "Coriolis" RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    # nu = size(RHSc,1)
    # @time "Viscous" if ν != 0
    #     RHSv = VB.rhs_viscosity(N,m; ns, ν)
    # else
    #     RHSv = spzeros(Complex{T}, nu, nu)
    # end

    RHSuu = RHSc #+ RHSv


    #wigner symbol temporary arrays alloc
    @everywhere begin
        Limace.DiscretePart.wig_table_init($(2N), 9)
        Limace.DiscretePart.wig_temp_init($(2N))
    end

    if B0poloidal
        @time "Lorentz" RHSub = DP.rhs_lorentz_bpol_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_bpol_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
        if mb0 != 0
            @time "Lorentz" RHSub += (-1)^mb0*DP.rhs_lorentz_bpol_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            @time "Induction" RHSbu += (-1)^mb0*DP.rhs_induction_bpol_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        @time "Lorentz" RHSub = DP.rhs_lorentz_btor_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
        @time "Induction" RHSbu = DP.rhs_induction_btor_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
        if mb0 != 0
            @time "Lorentz" RHSub += (-1)^mb0*DP.rhs_lorentz_btor_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            @time "Induction" RHSbu += (-1)^mb0*DP.rhs_induction_btor_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end

    #wigner symbol temporary arrays dealloc
    @everywhere Limace.DiscretePart.wig_temp_free()

    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

function rhs_combined(N,m; ns = 0, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal::Vector{Bool} = [true], lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], B0fac = [1.0], thresh = sqrt(eps()), kwargs...) where T
    
    @time "Coriolis" RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    # nu = size(RHSc,1)
    # @time "Viscous" if ν != 0
    #     RHSv = VB.rhs_viscosity(N,m; ns, ν)
    # else
    #     RHSv = spzeros(Complex{T}, nu, nu)
    # end

    RHSuu = RHSc #+ RHSv

    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)

    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

  
    nu = size(RHSuu,1)
    nb = size(RHSbb,1)

    RHSub = spzeros(ComplexF64,nu,nb)
    RHSbu = spzeros(ComplexF64,nb,nu)

    for (lmnb0,b0p,α) in zip(lmnb0,B0poloidal,B0fac)
        lb0,mb0,nb0 = lmnb0
        if b0p
            @time "Lorentz" RHSubt = DP.rhs_lorentz_bpol(N,m, lmnb0; ns, η, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_bpol(N,m, lmnb0; ns, η, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_bpol(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_bpol(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        else
            @time "Lorentz" RHSubt = DP.rhs_lorentz_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_btor(N,m, lmnb0; ns, η, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_btor(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_btor(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        end
        RHSub +=RHSubt*α
        RHSbu +=RHSbut*α
    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

   

    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

function rhs_dist_combined(N,m; ns = 0, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal::Vector{Bool} = [true], lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], B0fac = [1.0], thresh = sqrt(eps()), kwargs...) where T
    
    @time "Coriolis" RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    @time "Magnetic diffusion" RHSbb = BB.rhs_diffusion(N,m; ns, η)
    # nu = size(RHSc,1)
    # @time "Viscous" if ν != 0
    #     RHSv = VB.rhs_viscosity(N,m; ns, ν)
    # else
    #     RHSv = spzeros(Complex{T}, nu, nu)
    # end

    RHSuu = RHSc #+ RHSv


    #wigner symbol temporary arrays alloc
    @everywhere begin
        Limace.DiscretePart.wig_table_init($(2N), 9)
        Limace.DiscretePart.wig_temp_init($(2N))
    end

    nu = size(RHSuu,1)
    nb = size(RHSbb,1)
    RHSub = spzeros(ComplexF64,nu,nb)
    RHSbu = spzeros(ComplexF64,nb,nu)
    
    for (lmnb0,b0p,α) in zip(lmnb0,B0poloidal,B0fac)
        lb0,mb0,nb0 = lmnb0
        if b0p
            @time "Lorentz" RHSubt = DP.rhs_lorentz_bpol_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_bpol_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_bpol_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_bpol_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        else
            @time "Lorentz" RHSubt = DP.rhs_lorentz_btor_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_btor_dist(N,m, lmnb0; ns, η, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_btor_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_btor_dist(N,m, (lb0,-mb0,nb0); ns, η, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        end
        RHSub +=RHSubt*α
        RHSbu +=RHSbut*α
    end

    #wigner symbol temporary arrays dealloc
    @everywhere Limace.DiscretePart.wig_temp_free()


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
    

    lb0,mb0,nb0 = lmnb0
    RHSc = VB.rhs_coriolis(N,m; ns, Ω)

    RHSuu = RHSc #+ RHSv


    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

    if B0poloidal
        # RHSub = DP.rhs_lorentz_bpol_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
        # RHSbu = DP.rhs_induction_bpol_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
        # if mb0 != 0
        #     RHSub += (-1)^mb0*DP.rhs_lorentz_bpol_cond(N,m, (lb0,-mb0,nb0); ns, η=Ω, thresh, kwargs...)*B0fac
        #     RHSbu += (-1)^mb0*DP.rhs_induction_bpol_cond(N,m, (lb0,-mb0,nb0); ns, η=Ω, thresh, kwargs...)*B0fac
        # end
        # RHSub/=2
        # RHSbu/=2
    else
        RHSub = DP.rhs_lorentz_btor_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
        RHSbu = DP.rhs_induction_btor_cond(N,m, lmnb0; ns, η=Ω, thresh, kwargs...)*B0fac
        if mb0 != 0
            RHSub += DP.rhs_lorentz_btor_cond(N,m, (lb0,-mb0,nb0); ns, η=Ω, thresh, kwargs...)*B0fac*(-1)^mb0
            RHSbu += DP.rhs_induction_btor_cond(N,m, (lb0,-mb0,nb0); ns, η=Ω, thresh, kwargs...)*B0fac*(-1)^mb0
            RHSub/=2
            RHSbu/=2
        end

    end

    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()

    RHS = [RHSuu RHSub
           RHSbu spzeros(size(RHSbu,1),size(RHSub,2))]
    
    return RHS
end
end