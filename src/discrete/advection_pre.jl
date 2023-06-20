#poloidal B1, poloidal B2
const _advection_sss_pre = _lorentz_SSs_pre

#poloidal B1, toroidal B0
const _advection_sts_pre = _lorentz_STs_pre

#toroidal B1, toroidal B0
const _advection_tts_pre = _lorentz_TTs_pre

#poloidal B1, poloidal B0
const _advection_sst_pre = _lorentz_SSt_pre

#poloidal B1, toroidal B0
const _advection_stt_pre = _lorentz_STt_pre

#toroidal B1, toroidal B0
const _advection_ttt_pre = _lorentz_TTt_pre

#matrix assembly

function rhs_advection_upol_pre(N, m, lmnu0, r, wr, js_a1, js_a0;
    ns=false,
    η::T=1.0,
    thresh=sqrt(eps()),
    su0::Sf=s_in_pre,
    d_su0::dSf=d_s_in_pre,
    d2_su0::d2Sf=d2_s_in_pre,
    d3_su0::d3Sf=d3_s_in_pre,
    s_mf_u0::Sf2 = s_in,
    conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lu0, mu0, nu0 = lmnu0
    rls = [r .^ l for l in 1:N]

    Base.@propagate_inbounds Su0(l, m, n, r, i) = su0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSu0(l, m, n, r, i) = d_su0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Su0(l, m, n, r, i) = d2_su0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Su0(l, m, n, r, i) = d3_su0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds tu(l, m, n, r, i) = t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dtu(l, m, n, r, i) = d_t_in_pre(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2tu(l, m, n, r, i) = d2_t_in_pre(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds su(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dsu(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2su(l, m, n, r, i) = d2_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3su(l, m, n, r, i) = d3_s_in_pre(js_a1, rls, l, m, n, r, i)

    lmn_p = Limace.InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    is, js, aijs = Int[], Int[], Complex{T}[]


    @inbounds for (i, lmni) in enumerate(lmn_p)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && !ncondition(lu0, ni, nu0 + 1, nj) && continue
            conditions && !condition1(li, lu0, lj, mi, mu0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _advection_sss_pre(lmnj, lmnu0, lmni, r, wr, su, Su0, su, dsu, d2su, d3su, dSu0)
            aij += _advection_sss_pre(lmnu0, lmnj, lmni, r, wr, Su0, su, su, dSu0, d2Su0, d3Su0, dsu)
            appendit!(is, js, aijs, i, j, aij; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj, mj, nj = lmnj
            conditions && !ncondition(lu0, ni, nu0, nj) && continue
            conditions && !condition2(li, lu0, lj, mi, mu0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _advection_sts_pre(lmnu0, lmnj, lmni, r, wr, Su0, tu, su, dSu0, d2Su0, dtu, d2tu)
            appendit!(is, js, aijs, i, j + np, aij; thresh)

        end
    end

    @inbounds for (i, lmni) in enumerate(lmn_t)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && !ncondition(lu0, ni, nu0, nj) && continue
            conditions && !condition2(li, lu0, lj, mi, mu0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _advection_sst_pre(lmnu0, lmnj, lmni, r, wr, Su0, su, tu, dSu0, d2Su0)
            aij += _advection_sst_pre(lmnj, lmnu0, lmni, r, wr, su, Su0, tu, dsu, d2su)
            appendit!(is, js, aijs, i + np, j, aij; thresh)

        end
        for (j, lmnj) in enumerate(lmn_t)
            lj, mj, nj = lmnj
            conditions && !ncondition(lu0, ni, nu0, nj) && continue
            conditions && !condition1(li, lu0, lj, mi, mu0, mj) && continue
            # _dummy!(is,js,aijs,i,j)
            aij = _advection_stt_pre(lmnu0, lmnj, lmni, r, wr, Su0, tu, tu, dSu0, dtu)
            appendit!(is, js, aijs, i + np, j + np, aij; thresh)

        end
    end
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatu, nmatu)
end

function rhs_advection_utor_pre(N, m, lmnu0, r, wr, js_a1, js_a0;
    ns=false,
    η::T=1.0,
    thresh=sqrt(eps()),
    tu0::Tf=t_in_pre,
    d_tu0::dTf=d_t_in_pre,
    d2_tu0::d2Tf=d2_t_in_pre,
    d3_tu0::d3Tf=d3_t_in_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    rhs_lorentz_btor_cond_pre(N, m, lmnu0, r, wr, js_a1, js_a0; conditions, ns, η, thresh,
    tmfb0=tu0,
    d_tmfb0=d_tu0,
    d2_tmfb0=d2_tu0,
    d3_tmfb0=d3_tu0)
end
