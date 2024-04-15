module MHDProblem

using Limace
using SparseArrays
using Wigxjpf
using Distributed

const VB = Limace.InviscidBasis
const BB = Limace.InsulatingMFBasis
const DP = Limace.DiscretePart 

function lhs(N,m; ns = false, Ω::T = 1.0, η::T = 1.0, external=true) where T
    LHSu = VB.lhs(N,m; ns, Ω)
    if external
        LHSb = BB.lhs(N,m; ns, Ω)
    else
        LHSb = DP.lhs_inertial_b(N,m)
    end

    LHS = blockdiag(LHSu,LHSb)
    return LHS
end

function rhs(N, m; 
                ns = false, 
                Ω::T = 2.0, 
                ν::T = 1.0, 
                η::T = 1.0, 
                B0poloidal = true, 
                lmnb0::NTuple{3,Int} = (1,0,1), 
                B0fac = 1.0,
                thresh = 1000eps(), 
                u0poloidal = true,
                u0fac = 0.0,
                ufunc = DP.t_in,
                external = true,
                lmnu0::NTuple{3,Int} = (1,0,1),
                kwargs...) where T
    
    lb0,mb0,nb0 = lmnb0
    RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    if external
        RHSbb = BB.rhs_diffusion(N,m; ns, η)
    else
        RHSbb = DP.rhs_diffusion(N,m; η)
    end
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
        RHSub = DP.rhs_lorentz_bpol(N,m, lmnb0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_bpol(N,m, lmnb0; ns, thresh, external, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_bpol(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_bpol(N,m, (lb0,-mb0,nb0); ns, thresh, external, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        RHSub = DP.rhs_lorentz_btor(N,m, lmnb0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_btor(N,m, lmnb0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_btor(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_btor(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end

    RHSub*=B0fac
    RHSbu*=B0fac
    
    if u0fac != 0.0
        lu0, mu0, nu0 = lmnu0
        if u0poloidal
            RHSbbt = DP.rhs_induction_upol(N,m, lmnu0; ns, thresh, su = ufunc, kwargs...)
            if mu0 != 0
                RHSbbt += (-1)^mu0*DP.rhs_induction_upol(N,m, (lu0,-mu0,nu0); ns, thresh, su = ufunc, kwargs...)
                RHSbbt/=2
            end
        else
            RHSbbt = DP.rhs_induction_utor(N,m, lmnu0; ns, thresh, tu = ufunc, kwargs...)
            if mu0 != 0
                RHSbbt += (-1)^mu0*DP.rhs_induction_utor(N,m, (lu0,-mu0,nu0); ns, thresh, tu = ufunc, kwargs...)
                RHSbbt/=2
            end
        end
        RHSbb += u0fac*RHSbbt
    end



    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end


function rhs_pre(N, m; 
    ns = false, 
    Ω::T = 2.0, 
    ν::T = 1.0, 
    η::T = 1.0, 
    B0poloidal = true, 
    lmnb0::NTuple{3,Int} = (1,0,1), 
    B0fac=1.0,
    thresh = 1000eps(), 
    u0poloidal = true,
    u0fac = 0.0,
    ufunc = DP.t_in,
    nr_extra = 5,
    lmnu0::NTuple{3,Int} = (1,0,1),
    kwargs...) where T

    lb0,mb0,nb0 = lmnb0
    r,wr = DP.rquad(N+lb0+nb0+nr_extra)

	js_a1 = DP.jacobis_l(N,r,1.0)
	js_a0 = DP.jacobis_l(N,r,0.0)


    RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    RHSbb = BB.rhs_diffusion(N,m; ns, η)
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
        RHSub = DP.rhs_lorentz_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_bpol_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_bpol_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        RHSub = DP.rhs_lorentz_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_btor_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_btor_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end
    RHSbu*=B0fac
    RHSub*=B0fac

    # if u0fac != 0.0
    #     lu0, mu0, nu0 = lmnu0
    #     if u0poloidal
    #         RHSbbt = DP.rhs_induction_upol(N,m, lmnu0; ns, thresh, su = ufunc, kwargs...)
    #         if mu0 != 0
    #             RHSbbt += (-1)^mu0*DP.rhs_induction_upol(N,m, (lu0,-mu0,nu0); ns, thresh, su = ufunc, kwargs...)
    #             RHSbbt/=2
    #         end
    #     else
    #         RHSbbt = DP.rhs_induction_utor(N,m, lmnu0; ns, thresh, tu = ufunc, kwargs...)
    #         if mu0 != 0
    #             RHSbbt += (-1)^mu0*DP.rhs_induction_utor(N,m, (lu0,-mu0,nu0); ns, thresh, tu = ufunc, kwargs...)
    #             RHSbbt/=2
    #         end
    #     end
    #     RHSbb += u0fac*RHSbbt
    # end



    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]

    return RHS
end

function rhs_pre_dist(N, m; 
    ns = false, 
    Ω::T = 2.0, 
    ν::T = 1.0, 
    η::T = 1.0, 
    B0poloidal = true, 
    lmnb0::NTuple{3,Int} = (1,0,1), 
    thresh = 1000eps(), 
    u0poloidal = true,
    u0fac = 0.0,
    ufunc = DP.t_in,
    lmnu0::NTuple{3,Int} = (1,0,1),
    kwargs...) where T

    lb0,mb0,nb0 = lmnb0
    r,wr = DP.rquad(N+lb0+nb0+5)

	@time js_a1 = pmap(l->DP.jacobis(N,1.0,l+1/2,r),1:N)
	@time js_a0 = pmap(l->DP.jacobis(N,0.0,l+1/2,r), 1:N)


    @time RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    @time RHSbb = BB.rhs_diffusion(N,m; ns, η)
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
        @time RHSub = DP.rhs_lorentz_bpol_dist_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        @time RHSbu = DP.rhs_induction_bpol_dist_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_bpol_dist_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_bpol_dist_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        RHSub = DP.rhs_lorentz_btor_dist_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_btor_dist_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_btor_dist_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_btor_dist_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end

    # if u0fac != 0.0
    #     lu0, mu0, nu0 = lmnu0
    #     if u0poloidal
    #         RHSbbt = DP.rhs_induction_upol(N,m, lmnu0; ns, thresh, su = ufunc, kwargs...)
    #         if mu0 != 0
    #             RHSbbt += (-1)^mu0*DP.rhs_induction_upol(N,m, (lu0,-mu0,nu0); ns, thresh, su = ufunc, kwargs...)
    #             RHSbbt/=2
    #         end
    #     else
    #         RHSbbt = DP.rhs_induction_utor(N,m, lmnu0; ns, thresh, tu = ufunc, kwargs...)
    #         if mu0 != 0
    #             RHSbbt += (-1)^mu0*DP.rhs_induction_utor(N,m, (lu0,-mu0,nu0); ns, thresh, tu = ufunc, kwargs...)
    #             RHSbbt/=2
    #         end
    #     end
    #     RHSbb += u0fac*RHSbbt
    # end



    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]

    return RHS
end


function rhs_dist(N,m; ns = false, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal = true, lmnb0::NTuple{3,Int} = (1,0,1), thresh = 1000eps(), kwargs...) where T
    
    lb0,mb0,nb0 = lmnb0
    RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    RHSbb = BB.rhs_diffusion(N,m; ns, η)
    
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
        RHSub = DP.rhs_lorentz_bpol_dist(N,m, lmnb0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_bpol_dist(N,m, lmnb0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_bpol_dist(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_bpol_dist(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    else
        RHSub = DP.rhs_lorentz_btor_dist(N,m, lmnb0; ns, thresh, kwargs...)
        RHSbu = DP.rhs_induction_btor_dist(N,m, lmnb0; ns, thresh, kwargs...)
        if mb0 != 0
            RHSub += (-1)^mb0*DP.rhs_lorentz_btor_dist(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSbu += (-1)^mb0*DP.rhs_induction_btor_dist(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
            RHSub/=2
            RHSbu/=2
        end
    end

    #wigner symbol temporary arrays dealloc
    @everywhere Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end

function rhs_combined(N,m; ns = false, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal::Vector{Bool} = [true], lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], B0fac = [1.0], thresh = sqrt(eps()), kwargs...) where T
    
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
            @time "Lorentz" RHSubt = DP.rhs_lorentz_bpol(N,m, lmnb0; ns, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_bpol(N,m, lmnb0; ns, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_bpol(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_bpol(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        else
            @time "Lorentz" RHSubt = DP.rhs_lorentz_btor(N,m, lmnb0; ns, thresh, kwargs...)
            @time "Induction" RHSbut = DP.rhs_induction_btor(N,m, lmnb0; ns, thresh, kwargs...)
            if mb0 != 0
                @time "Lorentz" RHSubt += (-1)^mb0*DP.rhs_lorentz_btor(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
                @time "Induction" RHSbut += (-1)^mb0*DP.rhs_induction_btor(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
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

function rhs_pre_combined(N, m; 
    ns = false, 
    Ω::T = 2.0, 
    η::T = 1.0, 
    B0poloidal::Vector{Bool} = [true], 
    lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], 
    B0fac = [1.0],
    thresh = 1000eps(),
    kwargs...) where T

    lb0,mb0,nb0 = [maximum(getindex.(lmnb0,i)) for i=1:3]
    r,wr = DP.rquad(N+lb0+nb0+5)

	js_a1 = DP.jacobis_l(N,r,1.0)
	js_a0 = DP.jacobis_l(N,r,0.0)


    RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHSuu = RHSc


    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

    nu = size(RHSuu,1)
    nb = size(RHSbb,1)
    
    RHSub = spzeros(ComplexF64,nu,nb)
    RHSbu = spzeros(ComplexF64,nb,nu)

    for (lmnb0,B0poloidal,B0fac) in zip(lmnb0,B0poloidal,B0fac)
        lb0,mb0,nb0 = lmnb0
        if B0poloidal
            RHSubt = DP.rhs_lorentz_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbut = DP.rhs_induction_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            if mb0 != 0
                RHSubt += (-1)^mb0*DP.rhs_lorentz_bpol_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
                RHSbut += (-1)^mb0*DP.rhs_induction_bpol_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        else
            RHSubt = DP.rhs_lorentz_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbut = DP.rhs_induction_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            if mb0 != 0
                RHSubt += (-1)^mb0*DP.rhs_lorentz_btor_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
                RHSbut += (-1)^mb0*DP.rhs_induction_btor_pre(N,m, (lb0,-mb0,nb0),r,wr, js_a1,js_a0; ns, thresh, kwargs...)
                RHSubt/=2
                RHSbut/=2
            end
        end
        RHSbu += RHSbut*B0fac
        RHSub += RHSubt*B0fac
    end



    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]

    return RHS
end

function rhs_pre_combined_noextra(N, m; 
    ns = false, 
    Ω::T = 2.0, 
    η::T = 1.0, 
    B0poloidal::Vector{Bool} = [true], 
    lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], 
    B0fac = [1.0],
    thresh = 1000eps(),
    kwargs...) where T

    lb0,mb0,nb0 = [maximum(getindex.(lmnb0,i)) for i=1:3]
    r,wr = DP.rquad(N+lb0+nb0+5)

	js_a1 = DP.jacobis_l(N,r,1.0)
	js_a0 = DP.jacobis_l(N,r,0.0)


    RHSc = VB.rhs_coriolis(N,m; ns, Ω)
    RHSbb = BB.rhs_diffusion(N,m; ns, η)

    RHSuu = RHSc


    #wigner symbol temporary arrays alloc
    wig_table_init(2N, 9)
    wig_temp_init(2N)

    nu = size(RHSuu,1)
    nb = size(RHSbb,1)
    
    RHSub = spzeros(ComplexF64,nu,nb)
    RHSbu = spzeros(ComplexF64,nb,nu)

    for (lmnb0,B0poloidal,B0fac) in zip(lmnb0,B0poloidal,B0fac)
        lb0,mb0,nb0 = lmnb0
        if B0poloidal
            RHSubt = DP.rhs_lorentz_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbut = DP.rhs_induction_bpol_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        else
            RHSubt = DP.rhs_lorentz_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
            RHSbut = DP.rhs_induction_btor_pre(N,m, lmnb0,r,wr, js_a1,js_a0; ns, thresh, kwargs...)
        end
        RHSbu += RHSbut*B0fac
        RHSub += RHSubt*B0fac
    end



    #wigner symbol temporary arrays dealloc
    Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]

    return RHS
end

function rhs_lorentz_induction(N,m,b0p,lmnb0; ns, thresh, kwargs...)
    lb0,mb0,nb0 = lmnb0
    florentz = b0p ? DP.rhs_lorentz_bpol_dist : DP.rhs_lorentz_btor_dist
    finduction = b0p ? DP.rhs_induction_bpol_dist : DP.rhs_induction_btor_dist
    @time "Lorentz $lb0, $mb0, $nb0" RHSubt = florentz(N,m, lmnb0; ns, thresh, kwargs...)
    @time "Induction $lb0, $mb0, $nb0" RHSbut = finduction(N,m, lmnb0; ns, thresh, kwargs...)
    if mb0 != 0
        @time "Lorentz $lb0, -$mb0, $nb0" RHSubt2 = florentz(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
        @time "Induction $lb0, -$mb0, $nb0" RHSbut2 = finduction(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)
        if mb0<0
            RHSubt = (RHSubt - (-1)^mb0*RHSubt2)*im/sqrt(2)
            RHSbut = (RHSbut - (-1)^mb0*RHSbut2)*im/sqrt(2)
        else
            RHSubt = (RHSubt2 + (-1)^mb0*RHSubt)/sqrt(2)
            RHSbut = (RHSbut2 + (-1)^mb0*RHSbut)/sqrt(2)
        end
    end
    return RHSubt, RHSbut
end

function rhs_dist_combined(N,m; ns = false, Ω::T = 2.0, ν::T = 1.0, η::T = 1.0, B0poloidal::Vector{Bool} = [true], lmnb0::Vector{NTuple{3,Int}} = [(1,0,1)], B0fac = [1.0], thresh = sqrt(eps()), kwargs...) where T
    
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
        RHSubt, RHSbut = rhs_lorentz_induction(N,m,b0p,lmnb0; ns, thresh, kwargs...)
        RHSub +=RHSubt*α
        RHSbu +=RHSbut*α
    end

    #wigner symbol temporary arrays dealloc
    @everywhere Limace.DiscretePart.wig_temp_free()


    RHS = [RHSuu RHSub
           RHSbu RHSbb]
    
    return RHS
end


function lhs_cond(N,m; ns = false, Ω::T = 1.0, η::T = 1.0) where T
    LHSu = VB.lhs(N,m; ns, Ω)
    LHSb = VB.lhs(N,m; ns, Ω)

    LHS = blockdiag(LHSu,LHSb)
    return LHS
end

function rhs_cond(N,m; ns = false, Ω::T = 2.0, B0poloidal = false, B0fac=1.0, lmnb0::NTuple{3,Int} = (1,0,1), thresh = sqrt(eps()), kwargs...) where T
    

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
        RHSub = DP.rhs_lorentz_btor_cond(N,m, lmnb0; ns, thresh, kwargs...)*B0fac
        RHSbu = DP.rhs_induction_btor_cond(N,m, lmnb0; ns, thresh, kwargs...)*B0fac
        if mb0 != 0
            RHSub += DP.rhs_lorentz_btor_cond(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)*B0fac*(-1)^mb0
            RHSbu += DP.rhs_induction_btor_cond(N,m, (lb0,-mb0,nb0); ns, thresh, kwargs...)*B0fac*(-1)^mb0
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

function rhs_cond_pre(N,m; ns = false, Ω::T = 2.0, B0poloidal = false, B0fac=1.0, lmnb0::NTuple{3,Int} = (1,0,1), thresh = sqrt(eps()), kwargs...) where T
    

    lb0,mb0,nb0 = lmnb0
    r,wr = DP.rquad(N+lb0+nb0+5)

	js_a1 = DP.jacobis_l(N,r,1.0)
	js_a0 = DP.jacobis_l(N,r,0.0)


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
        RHSub = DP.rhs_lorentz_btor_cond_pre(N,m, lmnb0, r, wr, js_a1, js_a0; ns, thresh, kwargs...)*B0fac
        RHSbu = DP.rhs_induction_btor_cond_pre(N,m, lmnb0, r, wr, js_a1, js_a0; ns, thresh, kwargs...)*B0fac
        if mb0 != 0
            RHSub += DP.rhs_lorentz_btor_cond_pre(N,m, (lb0,-mb0,nb0), r, wr, js_a1, js_a0; ns, thresh, kwargs...)*B0fac*(-1)^mb0
            RHSbu += DP.rhs_induction_btor_cond_pre(N,m, (lb0,-mb0,nb0), r, wr, js_a1, js_a0; ns, thresh, kwargs...)*B0fac*(-1)^mb0
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