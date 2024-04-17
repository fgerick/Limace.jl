##
## poloidal Lorentz force
##

Base.@propagate_inbounds function _D(S, dS, d2S, l, rgrid::Vector{T}, i) where {T}
    @inbounds r = rgrid[i]
    d2S(i) + 2 / r * dS(i) - l * (l + 1) / r^2 * S(i)
end
Base.@propagate_inbounds function _dD(S, dS, d2S, d3S, l, rgrid::Vector{T}, i) where {T}
    @inbounds r = rgrid[i]
    d3S(i) + 2 / r * d2S(i) - (2 + l * (l + 1)) / r^2 * dS(i) + 2l * (l + 1) / r^3 * S(i)
end
Base.@propagate_inbounds function _innert(t, t2, l, i)
    return l * (l + 1) * t(i) * t2(i)
end



#poloidal B1, poloidal B2
Base.@propagate_inbounds function _lorentz_SSs_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Sa, Sb, sc, dSa, d2Sa, d3Sa, dSb) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Aabc = adamgaunt(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Sa(i) = Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _dSa(i) = dSa(la, ma, na, r, i)
    Base.@propagate_inbounds _d2Sa(i) = d2Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _d3Sa(i) = d3Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _Sb(i) = Sb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _dSb(i) = dSb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _sc(i) = sc(lc, mc, nc, r, i)

    Base.@propagate_inbounds f1(i) = (p(lc) * (p(la) + p(lb) - p(lc)) * _D(_Sa, _dSa, _d2Sa, la, r, i) * (_Sb(i) + r[i] * _dSb(i)) +
                                      p(lb) * (p(la) - p(lb) + p(lc)) * r[i] * (_dD(_Sa, _dSa, _d2Sa, _d3Sa, la, r, i) * _Sb(i) + _D(_Sa, _dSa, _d2Sa, la, r, i) * _dSb(i))) / (2r[i]^2 * p(lc))

    Base.@propagate_inbounds f(i) = -_innert(_sc, f1, lc, i)
    # Base.@propagate_inbounds  f = r->_sc(i)*f1(i) #inners(_sc,f1,lc,r)

    aij = ∫dr_pre(f, r, wr) * Aabc
    return aij
end


#poloidal B1, toroidal B0
Base.@propagate_inbounds function _lorentz_STs_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Sa, Tb, sc, dSa, d2Sa, dTb, d2Tb) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Eabc = elsasser(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Sa(i) = Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _dSa(i) = dSa(la, ma, na, r, i)
    Base.@propagate_inbounds _d2Sa(i) = d2Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _Tb(i) = Tb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _dTb(i) = dTb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _d2Tb(i) = d2Tb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _sc(i) = sc(lc, mc, nc, r, i)


    Base.@propagate_inbounds f1(i) = (p(lc) * r[i]^2 * _D(_Sa, _dSa, _d2Sa, la, r, i) * _Tb(i) +
                                      (p(la) + p(lb) + p(lc)) * _Sa(i) * _Tb(i) -
                                      (p(la) + p(lb) - p(lc)) * (r[i] * _Sa(i) * _dTb(i) +
                                                                 r[i] * _dSa(i) * _Tb(i) +
                                                                 r[i]^2 * _dSa(i) * _dTb(i)) -
                                      p(lb) * r[i]^2 * _d2Sa(i) * _Tb(i) -
                                      p(la) * r[i]^2 * _d2Tb(i) * _Sa(i)
    ) / (r[i]^3 * p(lc))

    Base.@propagate_inbounds f(i) = -innert(_sc, f1, lc, i)

    aij = ∫dr_pre(f, r, wr) * Eabc
    return aij
end


#toroidal B1, toroidal B0
Base.@propagate_inbounds function _lorentz_TTs_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Ta, Tb, sc, dTa, dTb) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Aabc = adamgaunt(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Ta(i) = Ta(la, ma, na, r, i)
    Base.@propagate_inbounds _dTa(i) = dTa(la, ma, na, r, i)
    Base.@propagate_inbounds _Tb(i) = Tb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _dTb(i) = dTb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _sc(i) = sc(lc, mc, nc, r, i)


    Base.@propagate_inbounds f1(i) = (p(lc) * (p(la) + p(lb) - p(lc)) * (_Ta(i) + r[i] * _dTa(i)) * _Tb(i) + p(la) * (-p(la) + p(lb) + p(lc)) * r[i] * (_Ta(i) * _dTb(i) + _dTa(i) * _Tb(i))) / (2r[i]^2 * p(lc))

    Base.@propagate_inbounds f(i) = -innert(_sc, f1, lc, i)

    aij = ∫dr_pre(f, r, wr) * Aabc
    return aij
end



##
## toroidal lorentz equation
##


#poloidal B1, poloidal B0
Base.@propagate_inbounds function _lorentz_SSt_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Sa, Sb, tc, dSa, d2Sa) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Eabc = elsasser(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Sa(i) = Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _dSa(i) = dSa(la, ma, na, r, i)
    Base.@propagate_inbounds _d2Sa(i) = d2Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _Sb(i) = Sb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _tc(i) = tc(lc, mc, nc, r, i)

    # Base.@propagate_inbounds  D(f,l,r) = ∂(r->∂(f,r),r) + 2/r * ∂(f,r) - l*(l+1)/r^2 *f(i)
    Base.@propagate_inbounds f1(i) = -p(lb) * _D(_Sa, _dSa, _d2Sa, la, r, i) * _Sb(i) / (r[i] * p(lc))

    Base.@propagate_inbounds f(i) = _innert(_tc, f1, lc, i)

    aij = ∫dr_pre(f, r, wr) * Eabc
    return aij
end

#poloidal B1, toroidal B0
Base.@propagate_inbounds function _lorentz_STt_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Sa, Tb, tc, dSa, dTb) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Aabc = adamgaunt(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Sa(i) = Sa(la, ma, na, r, i)
    Base.@propagate_inbounds _dSa(i) = dSa(la, ma, na, r, i)
    Base.@propagate_inbounds _Tb(i) = Tb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _dTb(i) = dTb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _tc(i) = tc(lc, mc, nc, r, i)


    Base.@propagate_inbounds f1(i) = (p(lb) * (p(lb) - p(la) - p(lc)) * (_Sa(i) + r[i] * _dSa(i)) * _Tb(i) - p(la) * (p(la) - p(lb) - p(lc)) * _Sa(i) * (_Tb(i) + r[i] * _dTb(i))) / (2r[i]^2 * p(lc))

    Base.@propagate_inbounds f(i) = innert(_tc, f1, lc, i)

    aij = ∫dr_pre(f, r, wr) * Aabc
    return aij
end


#toroidal B1, toroidal B0
Base.@propagate_inbounds function _lorentz_TTt_pre(lmna, lmnb, lmnc, r::Vector{T}, wr, Ta, Tb, tc) where {T}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    Eabc = elsasser(la, lb, lc, ma, mb, mc)
    Base.@propagate_inbounds _Ta(i) = Ta(la, ma, na, r, i)
    Base.@propagate_inbounds _Tb(i) = Tb(lb, mb, nb, r, i)
    Base.@propagate_inbounds _tc(i) = tc(lc, mc, nc, r, i)


    Base.@propagate_inbounds f1(i) = p(la) * _Ta(i) * _Tb(i) / (r[i] * p(lc))

    Base.@propagate_inbounds f(i) = innert(_tc, f1, lc, i)

    aij = ∫dr_pre(f, r, wr) * Eabc
    return aij
end

#create mutating functions for easier sparse matrices assembly
# flist = [:_lorentz_SSs_pre, :_lorentz_STs_pre, :_lorentz_TTs_pre,
#          :_lorentz_SSt_pre, :_lorentz_STt_pre, :_lorentz_TTt_pre]

# for f in flist
#     @eval begin
#         function $(Symbol(string(f)*"!"))(is::Vector{Int},js::Vector{Int},aijs::Vector{Complex{T}},i,j, lmna, lmnb, lmnc, r::Vector{T},wr, fa,fb,fc, args...; thresh=sqrt(eps())) where T
#             aij = $(f)(lmna, lmnb, lmnc, r,wr, fa,fb,fc, args...)
#             if abs(aij) > thresh
#                 push!(is,i)
#                 push!(js,j)
#                 push!(aijs,aij) 
#             end
#             return nothing
#         end
#     end
# end


#matrix assembly

function rhs_lorentz_bpol_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    smfb0::Sf=s_mf_pre,
    d_smfb0::dSf=d_s_mf_pre,
    d2_smfb0::d2Sf=d2_s_mf_pre,
    d3_smfb0::d3Sf=d3_s_mf_pre,
    s_mf_b0::Sf2 = s_mf,
    conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Smfb0(l, m, n, r, i) = smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmfb0(l, m, n, r, i) = d_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smfb0(l, m, n, r, i) = d2_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smfb0(l, m, n, r, i) = d3_smfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]


    # rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
    @inbounds for (i, lmni) in enumerate(lmn_p)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_SSs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, su, dSmf, d2Smf, d3Smf, dSmfb0)
            aij += _lorentz_SSs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, su, dSmfb0, d2Smfb0, d3Smfb0, dSmf)
            appendit!(is, js, aijs, i, j, aij; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_STs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, su, dSmfb0, d2Smfb0, dTmf, d2Tmf)
            appendit!(is, js, aijs, i, j + npb, aij; thresh)

        end
    end

    @inbounds for (i, lmni) in enumerate(lmn_t)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_SSt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, tu, dSmfb0, d2Smfb0)
            aij += _lorentz_SSt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, tu, dSmf, d2Smf)
            appendit!(is, js, aijs, i + np, j, aij; thresh)

        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_STt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, tu, dSmfb0, dTmf)
            appendit!(is, js, aijs, i + np, j + npb, aij; thresh)

        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function rhs_lorentz_bpol_pre_thread(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    smfb0::Sf=s_mf_pre,
    d_smfb0::dSf=d_s_mf_pre,
    d2_smfb0::d2Sf=d2_s_mf_pre,
    d3_smfb0::d3Sf=d3_s_mf_pre,
    s_mf_b0::Sf2 = s_mf,
    conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Smfb0(l, m, n, r, i) = smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmfb0(l, m, n, r, i) = d_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smfb0(l, m, n, r, i) = d2_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smfb0(l, m, n, r, i) = d3_smfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    nthreads = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:nthreads],[Int[] for _ in 1:nthreads], [Complex{T}[] for _ in 1:nthreads]


    # rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
    @inbounds Threads.@threads for i in eachindex(lmn_p)
        lmni = lmn_p[i]
        it = Threads.threadid()
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_SSs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, su, dSmf, d2Smf, d3Smf, dSmfb0)
            aij += _lorentz_SSs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, su, dSmfb0, d2Smfb0, d3Smfb0, dSmf)
            appendit!(is[it], js[it], aijs[it], i, j, aij; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_STs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, su, dSmfb0, d2Smfb0, dTmf, d2Tmf)
            appendit!(is[it], js[it], aijs[it], i, j + npb, aij; thresh)

        end
    end

    @inbounds Threads.@threads for i in eachindex(lmn_t)
        lmni = lmn_t[i]
        it = Threads.threadid()
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_SSt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, tu, dSmfb0, d2Smfb0)
            aij += _lorentz_SSt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, tu, dSmf, d2Smf)
            appendit!(is[it], js[it], aijs[it], i + np, j, aij; thresh)

        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _lorentz_STt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, tu, dSmfb0, dTmf)
            appendit!(is[it], js[it], aijs[it], i + np, j + npb, aij; thresh)

        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function rhs_lorentz_btor_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    tmfb0::Tf=t_mf_pre,
    d_tmfb0::dTf=d_t_mf_pre,
    d2_tmfb0::d2Tf=d2_t_mf_pre,
    d3_tmfb0::d3Tf=d3_t_mf_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Tmfb0(l, m, n, r, i) = tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmfb0(l, m, n, r, i) = d_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmfb0(l, m, n, r, i) = d2_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmfb0(l, m, n, r, i) = d3_tmfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]


    for (i, lmni) in enumerate(lmn_p)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, su, dSmf, d2Smf, dTmfb0, d2Tmfb0)
            appendit!(is, js, aijs, i, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTs_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, su, dTmf, dTmfb0)
            aij += _lorentz_TTs_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, su, dTmfb0, dTmf)
            appendit!(is, js, aijs, i, j + npb, aij; thresh)
        end
    end

    for (i, lmni) in enumerate(lmn_t)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, tu, dSmf, dTmfb0)
            appendit!(is, js, aijs, i + np, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTt_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, tu)
            aij += _lorentz_TTt_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, tu)
            appendit!(is, js, aijs, i + np, j + npb, aij; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function rhs_lorentz_btor_cond_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    tmfb0::Tf=t_in_pre,
    d_tmfb0::dTf=d_t_in_pre,
    d2_tmfb0::d2Tf=d2_t_in_pre,
    d3_tmfb0::d3Tf=d3_t_in_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Tmfb0(l, m, n, r, i) = tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmfb0(l, m, n, r, i) = d_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmfb0(l, m, n, r, i) = d2_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmfb0(l, m, n, r, i) = d3_tmfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_in_pre(js_a1, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InviscidBasis.lmn_upol(N, m, ns)
    lmn_bt = InviscidBasis.lmn_utor(N, m, ns)

    npb = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]


    for (i, lmni) in enumerate(lmn_p)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, su, dSmf, d2Smf, dTmfb0, d2Tmfb0)
            appendit!(is, js, aijs, i, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTs_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, su, dTmf, dTmfb0)
            aij += _lorentz_TTs_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, su, dTmfb0, dTmf)
            appendit!(is, js, aijs, i, j + npb, aij; thresh)
        end
    end

    for (i, lmni) in enumerate(lmn_t)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, tu, dSmf, dTmfb0)
            appendit!(is, js, aijs, i + np, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTt_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, tu)
            aij += _lorentz_TTt_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, tu)
            appendit!(is, js, aijs, i + np, j + npb, aij; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function rhs_lorentz_bpol_dist_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    smfb0::Sf=s_mf_pre,
    d_smfb0::dSf=d_s_mf_pre,
    d2_smfb0::d2Sf=d2_s_mf_pre,
    d3_smfb0::d3Sf=d3_s_mf_pre,
    s_mf_b0::Sf2 = s_mf,
    conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Smfb0(l, m, n, r, i) = smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmfb0(l, m, n, r, i) = d_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smfb0(l, m, n, r, i) = d2_smfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smfb0(l, m, n, r, i) = d3_smfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    nt = nprocs()
    isd,jsd,aijsd = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


    # rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
    @sync @distributed for i in shuffle(eachindex(lmn_p))
        lmni = lmn_p[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            # r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            aij = _lorentz_SSs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, su, dSmf, d2Smf, d3Smf, dSmfb0)
            aij += _lorentz_SSs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, su, dSmfb0, d2Smfb0, d3Smfb0, dSmf)
            appendit!(is, js, aijs, i, j, aij; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            # r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            aij = _lorentz_STs_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, su, dSmfb0, d2Smfb0, dTmf, d2Tmf)
            appendit!(is, js, aijs, i, j + npb, aij; thresh)

        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_t))
        lmni = lmn_t[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            # r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            aij = _lorentz_SSt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Smf, tu, dSmfb0, d2Smfb0)
            aij += _lorentz_SSt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Smfb0, tu, dSmf, d2Smf)
            appendit!(is, js, aijs, i + np, j, aij; thresh)

        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            # r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            aij = _lorentz_STt_pre(lmnb0, lmnj, lmni, r, wr, Smfb0, Tmf, tu, dSmfb0, dTmf)
            appendit!(is, js, aijs, i + np, j + npb, aij; thresh)

        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(vcat(isd...),vcat(jsd...),vcat(aijsd...), nmatu, nmatb)
end

function rhs_lorentz_btor_dist_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    thresh::T=sqrt(eps()),
    tmfb0::Tf=t_mf_pre,
    d_tmfb0::dTf=d_t_mf_pre,
    d2_tmfb0::d2Tf=d2_t_mf_pre,
    d3_tmfb0::d3Tf=d3_t_mf_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    # r, wr = rquad(N+lb0+nb0+5)
    rls = [r .^ l for l in 1:N]
    # js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
    # js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

    Base.@propagate_inbounds Tmfb0(l, m, n, r, i) = tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmfb0(l, m, n, r, i) = d_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmfb0(l, m, n, r, i) = d2_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmfb0(l, m, n, r, i) = d3_tmfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Smf(l, m, n, r, i) = d3_s_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Tmf(l, m, n, r, i) = t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmf(l, m, n, r, i) = d_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmf(l, m, n, r, i) = d2_t_mf_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmf(l, m, n, r, i) = d3_t_mf_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    nt = nprocs()
    isd,jsd,aijsd = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


    @sync @distributed for i in shuffle(eachindex(lmn_p))
        lmni = lmn_p[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STs_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, su, dSmf, d2Smf, dTmfb0, d2Tmfb0)
            appendit!(is, js, aijs, i, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTs_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, su, dTmf, dTmfb0)
            aij += _lorentz_TTs_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, su, dTmfb0, dTmf)
            appendit!(is, js, aijs, i, j + npb, aij; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_t))
        lmni = lmn_t[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_STt_pre(lmnj, lmnb0, lmni, r, wr, Smf, Tmfb0, tu, dSmf, dTmfb0)
            appendit!(is, js, aijs, i + np, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0+1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            aij = _lorentz_TTt_pre(lmnj, lmnb0, lmni, r, wr, Tmf, Tmfb0, tu)
            aij += _lorentz_TTt_pre(lmnb0, lmnj, lmni, r, wr, Tmfb0, Tmf, tu)
            appendit!(is, js, aijs, i + np, j + npb, aij; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(vcat(isd...),vcat(jsd...),vcat(aijsd...), nmatu, nmatb)
end