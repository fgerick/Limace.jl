# @inline p(l) = l*(l+1.0)

# _D(S,dS,d2S,l,rgrid::Vector{T},i)

##
## poloidal induction equation
##

# @inline function _inners(s,s2, l, r) 
#     return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r[i]*s(r),r)*∂(r->r[i]*s2(r),r))/r[i]^2
# end

@inline function _inners(s,s2, ds, ds2, l, r,i) 
    return l*(l+1)*(s(i)*s2(i)*l*(l+1)+(s(i)+r[i]*ds(i))*(s2(i)+r[i]*ds2(i)))/r[i]^2
end

# @inline function _inners(s,s2, ds2, d2s2, l, rgrid,i) 
#     return l*(l+1)*(s(i)*_D(s2,ds2,d2s2,l,rgrid, i))
# end

#poloidal flow, poloidal B0
function _induction_sSS_pre(lmna, lmnb, lmnc, r,wr, sa,Sb,Sc, dsa, dSb, dSc, d2sa, d2Sb, __sa, __Sb, __Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    Base.@propagate_inbounds _sa(i) = sa(la,ma,na,r,i)
    Base.@propagate_inbounds _dsa(i) = dsa(la,ma,na,r,i)
    Base.@propagate_inbounds _d2sa(i) = d2sa(la,ma,na,r,i)
    Base.@propagate_inbounds _Sb(i) = Sb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _dSb(i) = dSb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _d2Sb(i) = d2Sb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _Sc(i) = Sc(lc,mc,nc,r,i)
    Base.@propagate_inbounds _dSc(i) = dSc(lc,mc,nc,r,i)

    Base.@propagate_inbounds f1(i) = (-p(la)*(-p(la)+p(lb)+p(lc))*_sa(i)*(_Sb(i)+r[i]*_dSb(i)) + 
                      p(lb)*( p(la)-p(lb)+p(lc))*_Sb(i)*(_sa(i)+r[i]*_dsa(i)))/(2r[i]^2*p(lc))
    Base.@propagate_inbounds df1(i) = -2/r[i]*f1(i) + (-p(la)*(-p(la)+p(lb)+p(lc))*(_dsa(i)*(_Sb(i)+r[i]*_dSb(i)) + _sa(i)*(2*_dSb(i)+r[i]*_d2Sb(i))) + 
        p(lb)*( p(la)-p(lb)+p(lc))*(_dSb(i)*(_sa(i)+r[i]*_dsa(i))+ _Sb(i)*(2*_dsa(i) + r[i]*_d2sa(i)) ) )/(2r[i]^2*p(lc))
     
    
    Base.@propagate_inbounds f(i) = _inners(f1, _Sc, df1, _dSc, lc, r, i)

    aij = ∫dr_pre(f,r,wr)*Aabc

    aij +=lc*p(lb)*( p(la)-p(lb)+p(lc))*__Sb(lb,mb,nb,1.0)*∂(r->__sa(la,ma,na,r),1.0)*__Sc(lc,mc,nc,1.0)/2*Aabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTS_pre(lmna, lmnb, lmnc, r,wr, sa,Tb,Sc,dsa, dTb, dSc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    Base.@propagate_inbounds _sa(i) = sa(la,ma,na,r,i)
    Base.@propagate_inbounds _dsa(i) = dsa(la,ma,na,r,i)
    Base.@propagate_inbounds _Tb(i) = Tb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _dTb(i) = dTb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _Sc(i) = Sc(lc,mc,nc,r,i)
    Base.@propagate_inbounds _dSc(i) = dSc(lc,mc,nc,r,i)

    Base.@propagate_inbounds f1(i) = p(la)*_sa(i)*_Tb(i)/(r[i]*p(lc))
    Base.@propagate_inbounds df1(i) = -f1(i)/r[i] + p(la)*(_dsa(i)*_Tb(i) + _sa(i)*_dTb(i) )/(r[i]*p(lc))

    Base.@propagate_inbounds f(i) = _inners(f1,_Sc, df1, _dSc, lc, r, i)
    
    aij = ∫dr_pre(f,r,wr)*Eabc
    return aij
end


#toroidal flow, poloidal B0
function _induction_tSS_pre(lmna, lmnb, lmnc, r,wr, ta,Sb,Sc, dta, dSb, dSc, __ta, __Sb, __Sc) #need to evaluate functions at r=1 (not in gauss grid)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    Base.@propagate_inbounds _ta(i) = ta(la,ma,na,r,i)
    Base.@propagate_inbounds _dta(i) = dta(la,ma,na,r,i)
    Base.@propagate_inbounds _Sb(i) = Sb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _dSb(i) = dSb(lb,mb,nb,r,i)
    Base.@propagate_inbounds _Sc(i) = Sc(lc,mc,nc,r,i)
    Base.@propagate_inbounds _dSc(i) = dSc(lc,mc,nc,r,i)


    Base.@propagate_inbounds f1(i) = p(lb)*_ta(i)*_Sb(i)/(r[i]*p(lc))
    Base.@propagate_inbounds df1(i) = -f1(i)/r[i] + p(lb)*(_dta(i)*_Sb(i) + _ta(i)*_dSb(i) )/(r[i]*p(lc))

    Base.@propagate_inbounds f(i) = _inners(f1,_Sc, df1, _dSc, lc, r, i)
    
    aij = ∫dr_pre(f,r,wr)*Eabc
    #add contribution from external ∫dV (1<r<∞), 
    #if toroidal velocity is not 0 at r=11 
    aij += lc*p(lb)*__ta(la,ma,na,1.0)*__Sb(lb,mb,nb,1.0)*__Sc(lc,mc,nc,1.0)*Eabc

    return aij
end


#toroidal flow, toroidal B0 
#always 0


##
## toroidal induction equation
##

# @inline function innert(t,t2, l, r) 
#     return l*(l+1)*t(r)*t2(r)
# end


#poloidal flow, poloidal B0
function _induction_sST_pre(lmna, lmnb, lmnc, r,wr, sa,Sb,Tc,dsa,dSb,d2sa,d2Sb)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _sa(i) = sa(la,ma,na,r,i)
    @inline _dsa(i) = dsa(la,ma,na,r,i)
    @inline _d2sa(i) = d2sa(la,ma,na,r,i)
    @inline _Sb(i) = Sb(lb,mb,nb,r,i)
    @inline _dSb(i) = dSb(lb,mb,nb,r,i)
    @inline _d2Sb(i) = d2Sb(lb,mb,nb,r,i)
    @inline _Tc(i) = Tc(lc,mc,nc,r,i)


    @inline f1(i) = ((p(la)+p(lb)+p(lc))*_sa(i)*_Sb(i) - 
                       (p(la)+p(lb)-p(lc))*(r[i]*_dsa(i)*_Sb(i) + r[i]*_sa(i)*_dSb(i) + r[i]^2*_dsa(i)*_dSb(i)) - 
                       p(la)*r[i]^2*_sa(i)*_d2Sb(i) - p(lb)*r[i]^2*_Sb(i)*_d2sa(i))/(r[i]^3*p(lc))
    
    @inline f(i) =  _innert(_Tc,f1, lc, i)

    aij = ∫dr_pre(f,r,wr)*Eabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTT_pre(lmna, lmnb, lmnc, r,wr, sa,Tb,Tc,dsa, dTb)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _sa(i) = sa(la,ma,na,r,i)
    @inline _dsa(i) = dsa(la,ma,na,r,i)
    @inline _Tb(i) = Tb(lb,mb,nb,r,i)
    @inline _dTb(i) = dTb(lb,mb,nb,r,i)
    @inline _Tc(i) = Tc(lc,mc,nc,r,i)


    @inline f1(i) = (-p(lc)*(p(la)+p(lb)-p(lc))*(_sa(i)*_Tb(i) + r[i]*_dsa(i)*_Tb(i)) +
                        p(la)*(p(la)-p(lb)-p(lc))*(r[i]*_dsa(i)*_Tb(i) + r[i]*_sa(i)*_dTb(i)))/(2r[i]^2*p(lc))

    @inline f(i) =  _innert(_Tc,f1, lc, i)

    aij = ∫dr_pre(f,r,wr)*Aabc
    return aij
end 

#toroidal flow, poloidal B0
function _induction_tST_pre(lmna, lmnb, lmnc, r,wr, ta,Sb,Tc,dta,dSb)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _ta(i) = ta(la,ma,na,r,i)
    @inline _dta(i) = dta(la,ma,na,r,i)
    @inline _Sb(i) = Sb(lb,mb,nb,r,i)
    @inline _dSb(i) = dSb(lb,mb,nb,r,i)
    @inline _Tc(i) = Tc(lc,mc,nc,r,i)


    @inline f1(i) = (p(lc)*(p(la)+p(lb)-p(lc))*(_ta(i)*_Sb(i) + r[i]*_ta(i)*_dSb(i)) -
                       p(lb)*(p(lb)-p(la)-p(lc))*(r[i]*_dta(i)*_Sb(i) + r[i]*_ta(i)*_dSb(i)))/(2r[i]^2*p(lc))

    @inline f(i) =  _innert(_Tc,f1, lc, i)

    aij = ∫dr_pre(f,r,wr)*Aabc

    # aij += f(1.0)*Aabc*lb*lc
    return aij
end 

#toroidal flow, toroidal B0
function _induction_tTT_pre(lmna, lmnb, lmnc, r,wr, ta,Tb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _ta(i) = ta(la,ma,na,r,i)
    @inline _Tb(i) = Tb(lb,mb,nb,r,i)
    @inline _Tc(i) = Tc(lc,mc,nc,r,i)


    @inline f1(i) = _ta(i)*_Tb(i)/r[i]

    @inline f(i) =  innert(_Tc,f1, lc, i)

    aij = ∫dr_pre(f,r,wr)*Eabc

    # aij += f(1.0)*Eabc
    return aij
end 


#matrix assembly

function rhs_induction_bpol_pre(N,m, lmnb0, r,wr, js_a1,js_a0; 
        ns = false, 
        η::T=1.0, 
        thresh = sqrt(eps()),
        smfb0::Sf = s_mf_pre,
        d_smfb0::dSf = d_s_mf_pre, 
        d2_smfb0::d2Sf = d2_s_mf_pre, 
        d3_smfb0::d3Sf = d3_s_mf_pre,
        s_mf_b0::Sf2 = s_mf,
        conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lb0,mb0,nb0 = lmnb0
    rls = [r.^l for l in 1:N]

    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    Base.@propagate_inbounds  Smfb0(l,m,n,r,i) =smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dSmfb0(l,m,n,r,i) =d_smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2Smfb0(l,m,n,r,i) =d2_smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d3Smfb0(l,m,n,r,i) =d3_smfb0(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  Smf(l,m,n,r,i) = s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dSmf(l,m,n,r,i) = d_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Smf(l,m,n,r,i) = d2_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Smf(l,m,n,r,i) = d3_s_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  Tmf(l,m,n,r,i) = t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dTmf(l,m,n,r,i) = d_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Tmf(l,m,n,r,i) = d2_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Tmf(l,m,n,r,i) = d3_t_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  tu(l,m,n,r,i) = t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dtu(l,m,n,r,i) = d_t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2tu(l,m,n,r,i) = d2_t_in_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  su(l,m,n,r,i) = s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dsu(l,m,n,r,i) = d_s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2su(l,m,n,r,i) = d2_s_in_pre(js_a1,rls,l,m,n,r,i)


    is,js,aijs = Int[],Int[],Complex{T}[]


    @inbounds for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sSS_pre(lmnj, lmnb0, lmni, r, wr, su, Smfb0, Smf, dsu, dSmfb0, dSmf, d2su, d2Smfb0, s_in, s_mf_b0, s_mf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tSS_pre(lmnj, lmnb0, lmni, r, wr, tu, Smfb0, Smf, dtu, dSmfb0, dSmf, t_in, s_mf_b0, s_mf) #remember here we have the surface term specially evaluated!
            appendit!(is,js,aijs,i,j+np,aij; thresh)
        end
    end

    @inbounds for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sST_pre(lmnj, lmnb0, lmni, r, wr, su, Smfb0, Tmf, dsu, dSmfb0, d2su, d2Smfb0)
            appendit!(is,js,aijs,i+npb,j,aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tST_pre(lmnj, lmnb0, lmni, r, wr, tu, Smfb0, Tmf, dtu, dSmfb0)
            appendit!(is,js,aijs,i+npb,j+np,aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_btor_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    η::T=1.0,
    thresh=sqrt(eps()),
    tmfb0::Tf=t_mf_pre,
    d_tmfb0::dTf=d_t_mf_pre,
    d2_tmfb0::d2Tf=d2_t_mf_pre,
    d3_tmfb0::d3Tf=d3_t_mf_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    rls = [r .^ l for l in 1:N]

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

    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]


    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTS_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Smf, dsu, dTmfb0, dSmf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTT_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Tmf, dsu, dTmfb0)
            appendit!(is,js,aijs,i+npb,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tTT_pre(lmnj,lmnb0,lmni, r, wr, tu, Tmfb0, Tmf)
            appendit!(is,js,aijs,i+npb,j+np,aij; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_upol_pre(N,m, lmnu0, r,wr, js_a1,js_a0; ns = false, η::T=1.0, thresh = sqrt(eps()), conditions = true) where T
    # su::Sf = s_in_pre,
    # d_su::dSf = d_s_in_pre, 
    # d2_su::d2Sf = d2_s_in_pre, 
    # d3_su::d3Sf = d3_s_in_pre,
    # s_mf_u::Sf2 = s_in, condition = true
    # ) where T

    lu0,mu0,nu0 = lmnu0

    rls = [r .^ l for l in 1:N]


	Base.@propagate_inbounds  Smf(l,m,n,r,i) = s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dSmf(l,m,n,r,i) = d_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Smf(l,m,n,r,i) = d2_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Smf(l,m,n,r,i) = d3_s_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  Tmf(l,m,n,r,i) = t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dTmf(l,m,n,r,i) = d_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Tmf(l,m,n,r,i) = d2_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Tmf(l,m,n,r,i) = d3_t_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  su(l,m,n,r,i) = s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dsu(l,m,n,r,i) = d_s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2su(l,m,n,r,i) = d2_s_in_pre(js_a1,rls,l,m,n,r,i)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    # rwrs = [rquad(n+lu0+nu0+1) for n in 1:N]

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_sSS_pre(lmnu0,lmnj,lmni, r, wr, su, Smf, Smf, dsu, dSmf, dSmf, d2su, d2Smf, s_in, s_mf, s_mf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            conditions && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_sTS_pre(lmnu0,lmnj,lmni, r, wr, su, Tmf, Smf, dsu, dTmf, dSmf)
            appendit!(is,js,aijs,i,j+np,aij; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_sST_pre(lmnu0,lmnj,lmni, r, wr, su, Smf, Tmf, dsu, dSmf, d2su, d2Smf)
            appendit!(is,js,aijs,i+np,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_sTT_pre(lmnu0,lmnj,lmni, r, wr, su, Tmf, Tmf, dsu, dTmf)
            appendit!(is,js,aijs,i+np,j+np,aij; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end

function rhs_induction_utor_pre(N,m, lmnu0, r,wr, js_a1,js_a0; ns = false, η::T=1.0, thresh = sqrt(eps()), conditions=true) where T

    lu0,mu0,nu0 = lmnu0

    rls = [r .^ l for l in 1:N]

	Base.@propagate_inbounds  Smf(l,m,n,r,i) = s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dSmf(l,m,n,r,i) = d_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Smf(l,m,n,r,i) = d2_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Smf(l,m,n,r,i) = d3_s_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  Tmf(l,m,n,r,i) = t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dTmf(l,m,n,r,i) = d_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2Tmf(l,m,n,r,i) = d2_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d3Tmf(l,m,n,r,i) = d3_t_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  tu(l,m,n,r,i) = t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dtu(l,m,n,r,i) = d_t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2tu(l,m,n,r,i) = d2_t_in_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds  su(l,m,n,r,i) = s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  dsu(l,m,n,r,i) = d_s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds  d2su(l,m,n,r,i) = d2_s_in_pre(js_a1,rls,l,m,n,r,i)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]


    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_tSS_pre(lmnu0,lmnj,lmni, r, wr, tu, Smf, Smf, dtu, dSmf, dSmf, t_in, s_mf, s_mf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_tST_pre(lmnu0,lmnj,lmni, r, wr, tu, Smf, Tmf, dtu, dSmf)
            appendit!(is,js,aijs,i+np,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            conditions && !ncondition(lu0,ni,nu0+1,nj) && continue
            conditions && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            aij = _induction_tTT_pre(lmnu0,lmnj,lmni, r, wr, tu, Tmf, Tmf)
            appendit!(is,js,aijs,i+np,j+np,aij; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end

function rhs_induction_btor_cond_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    η::T=1.0,
    thresh=sqrt(eps()),
    tmfb0::Tf=t_in_pre,
    d_tmfb0::dTf=d_t_in_pre,
    d2_tmfb0::d2Tf=d2_t_in_pre,
    d3_tmfb0::d3Tf=d3_t_in_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    rls = [r .^ l for l in 1:N]

    Base.@propagate_inbounds Tmfb0(l, m, n, r, i) = tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds dTmfb0(l, m, n, r, i) = d_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Tmfb0(l, m, n, r, i) = d2_tmfb0(js_a0, rls, l, m, n, r, i)
    Base.@propagate_inbounds d3Tmfb0(l, m, n, r, i) = d3_tmfb0(js_a0, rls, l, m, n, r, i)

    Base.@propagate_inbounds Smf(l, m, n, r, i) = s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds dSmf(l, m, n, r, i) = d_s_in_pre(js_a1, rls, l, m, n, r, i)
    Base.@propagate_inbounds d2Smf(l, m, n, r, i) = d2_s_in_pre(js1a1, rls, l, m, n, r, i)
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

    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_bt = Limace.InviscidBasis.lmn_utor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]


    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTS_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Smf, dsu, dTmfb0, dSmf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTT_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Tmf, dsu, dTmfb0)
            appendit!(is,js,aijs,i+npb,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tTT_pre(lmnj,lmnb0,lmni, r, wr, tu, Tmfb0, Tmf)
            appendit!(is,js,aijs,i+npb,j+np,aij; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end
# function rhs_induction_btor_cond_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps())) where T
#     su = s_in
#     tu = t_in
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     tmf = t_in
#     smf = s_in
#     lmn_bp = Limace.InviscidBasis.lmn_upol(N,m,ns)
#     lmn_bt = Limace.InviscidBasis.lmn_utor(N,m,ns)

#     npb = length(lmn_bp)

#     is,js,aijs = Int[],Int[],Complex{T}[]

#     r, wr = rquad(N+lb0+nb0+5)

#     for (i,lmni) in enumerate(lmn_bp)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_sTS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, su, tmf, smf; thresh)
#         end
#     end

#     for (i,lmni) in enumerate(lmn_bt)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_sTT!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, su, tmf, tmf; thresh)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_tTT!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, tmf, tmf; thresh)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(is,js,aijs,nmatb, nmatu)
# end

function rhs_induction_bpol_dist_pre(N,m, lmnb0, r,wr, js_a1,js_a0; 
    ns = false, 
    η::T=1.0, 
    thresh = sqrt(eps()),
    smfb0::Sf = s_mf_pre,
    d_smfb0::dSf = d_s_mf_pre, 
    d2_smfb0::d2Sf = d2_s_mf_pre, 
    d3_smfb0::d3Sf = d3_s_mf_pre,
    s_mf_b0::Sf2 = s_mf,
    conditions=true) where {T,Sf,dSf,d2Sf,d3Sf,Sf2}

    lb0,mb0,nb0 = lmnb0
    rls = [r.^l for l in 1:N]

    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    Base.@propagate_inbounds  Smfb0(l,m,n,r,i) =smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dSmfb0(l,m,n,r,i) =d_smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2Smfb0(l,m,n,r,i) =d2_smfb0(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d3Smfb0(l,m,n,r,i) =d3_smfb0(js_a0,rls,l,m,n,r,i)

    Base.@propagate_inbounds  Smf(l,m,n,r,i) = s_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dSmf(l,m,n,r,i) = d_s_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2Smf(l,m,n,r,i) = d2_s_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d3Smf(l,m,n,r,i) = d3_s_mf_pre(js_a0,rls,l,m,n,r,i)

    Base.@propagate_inbounds  Tmf(l,m,n,r,i) = t_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dTmf(l,m,n,r,i) = d_t_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2Tmf(l,m,n,r,i) = d2_t_mf_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d3Tmf(l,m,n,r,i) = d3_t_mf_pre(js_a0,rls,l,m,n,r,i)

    Base.@propagate_inbounds  tu(l,m,n,r,i) = t_in_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dtu(l,m,n,r,i) = d_t_in_pre(js_a0,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2tu(l,m,n,r,i) = d2_t_in_pre(js_a0,rls,l,m,n,r,i)

    Base.@propagate_inbounds  su(l,m,n,r,i) = s_in_pre(js_a1,rls,l,m,n,r,i)
    Base.@propagate_inbounds  dsu(l,m,n,r,i) = d_s_in_pre(js_a1,rls,l,m,n,r,i)
    Base.@propagate_inbounds  d2su(l,m,n,r,i) = d2_s_in_pre(js_a1,rls,l,m,n,r,i)



    nt = nprocs()
    isd,jsd,aijsd = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


    @sync @distributed for i in shuffle(eachindex(lmn_bp))
        lmni = lmn_bp[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sSS_pre(lmnj, lmnb0, lmni, r, wr, su, Smfb0, Smf, dsu, dSmfb0, dSmf, d2su, d2Smfb0, s_in, s_mf_b0, s_mf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tSS_pre(lmnj, lmnb0, lmni, r, wr, tu, Smfb0, Smf, dtu, dSmfb0, dSmf, t_in, s_mf_b0, s_mf) #remember here we have the surface term specially evaluated!
            appendit!(is,js,aijs,i,j+np,aij; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_bt))
        lmni = lmn_bt[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sST_pre(lmnj, lmnb0, lmni, r, wr, su, Smfb0, Tmf, dsu, dSmfb0, d2su, d2Smfb0)
            appendit!(is,js,aijs,i+npb,j,aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tST_pre(lmnj, lmnb0, lmni, r, wr, tu, Smfb0, Tmf, dtu, dSmfb0)
            appendit!(is,js,aijs,i+npb,j+np,aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(vcat(isd...),vcat(jsd...),vcat(aijsd...),nmatb, nmatu)
end


function rhs_induction_btor_dist_pre(N, m, lmnb0, r, wr, js_a1, js_a0;
    ns=false,
    η::T=1.0,
    thresh=sqrt(eps()),
    tmfb0::Tf=t_mf_pre,
    d_tmfb0::dTf=d_t_mf_pre,
    d2_tmfb0::d2Tf=d2_t_mf_pre,
    d3_tmfb0::d3Tf=d3_t_mf_pre,
    conditions=true) where {T,Tf,dTf,d2Tf,d3Tf}

    lb0, mb0, nb0 = lmnb0
    rls = [r .^ l for l in 1:N]

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

    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    nt = nprocs()
    isd,jsd,aijsd = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


    @sync @distributed for i in shuffle(eachindex(lmn_bp))
        lmni = lmn_bp[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTS_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Smf, dsu, dTmfb0, dSmf)
            appendit!(is,js,aijs,i,j,aij; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_bt))
        lmni = lmn_bt[i]
        li,mi,ni = lmni
        is,js, aijs = first(localpart(isd)),first(localpart(jsd)),first(localpart(aijsd))
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_sTT_pre(lmnj,lmnb0,lmni, r, wr, su, Tmfb0, Tmf, dsu, dTmfb0)
            appendit!(is,js,aijs,i+npb,j,aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0+1,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            aij = _induction_tTT_pre(lmnj,lmnb0,lmni, r, wr, tu, Tmfb0, Tmf)
            appendit!(is,js,aijs,i+npb,j+np,aij; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(vcat(isd...),vcat(jsd...),vcat(aijsd...),nmatb, nmatu)
end

# function rhs_induction_bpol_dist_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps()), smfb0::Sf = s_mf) where {T,Sf}
#     su = s_in
#     tu = t_in
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
#     lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

#     npb = length(lmn_bp)


#     nt = nprocs()
#     is,js,aijs = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])

#     # r, wr = rquad(N+lb0+nb0+5)
#     rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
#     @sync @distributed for i in shuffle(eachindex(lmn_bp))
#         lmni = lmn_bp[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             # _dummy!(is,js,aijs,i,j)
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,lmnj, lmnb0, lmni, r, wr, su, smfb0, s_mf; thresh)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             # _dummy!(is,js,aijs,i,j+np)
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_tSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j+np,lmnj, lmnb0, lmni, r, wr, tu, smfb0, s_mf; thresh)
#         end
#     end

#     @sync @distributed for i in shuffle(eachindex(lmn_bt))
#         lmni = lmn_bt[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j,lmnj, lmnb0, lmni, r, wr, su, smfb0, t_mf; thresh)
#             # _dummy!(is,js,aijs,i+npb,j)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_tST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j+np,lmnj, lmnb0, lmni, r, wr, tu, smfb0, t_mf; thresh)
#             # _dummy!(is,js,aijs,i+npb,j+np)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatb, nmatu)
# end

function rhs_induction_btor_dist_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps())) where T
    su = s_in
    tu = t_in
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    nt = nprocs()
    is,js,aijs = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


    # r, wr = rquad(N+lb0+nb0+5)
    rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]

    @sync @distributed for i in shuffle(eachindex(lmn_bp))
        lmni = lmn_bp[i]
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, s_mf; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_bt))
        lmni = lmn_bt[i]
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTT!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, t_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tTT!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, t_mf, t_mf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatb, nmatu)
end
