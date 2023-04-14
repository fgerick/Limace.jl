# @inline p(l) = l*(l+1.0)

# _D(S,dS,d2S,l,rgrid::Vector{T},i)

##
## poloidal induction equation
##

# @inline function _inners(s,s2, l, r) 
#     return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r*s(r),r)*∂(r->r*s2(r),r))/r^2
# end

@inline function _inners(s,s2, ds, ds2, l, r,i) 
    return l*(l+1)*(s(i)*s2(i)*l*(l+1)+(s(i)+r[i]*ds(i))*(s2(i)+r[i]*ds2(i)))/r[i]^2
end

# @inline function _inners(s,s2, ds2, d2s2, l, rgrid,i) 
#     return l*(l+1)*(s(i)*_D(s2,ds2,d2s2,l,rgrid, i))
# end

#poloidal flow, poloidal B0
function _induction_sSS_pre(lmna, lmnb, lmnc, r,wr, sa,Sb,Sc, dsa, dSb, dSc, d2sa, d2Sb)
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
function _induction_sST_pre(lmna, lmnb, lmnc, r,wr, sa,Sb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _sa(i) = sa(la,ma,na,r)
    @inline _Sb(i) = Sb(lb,mb,nb,r)
    @inline _Tc(i) = Tc(lc,mc,nc,r)


    @inline f1(i) = ((p(la)+p(lb)+p(lc))*_sa(r)*_Sb(r) - 
                       (p(la)+p(lb)-p(lc))*(r*∂(_sa,r)*_Sb(r) + r*_sa(r)*∂(_Sb,r) + r^2*∂(_sa,r)*∂(_Sb,r)) - p(la)*r^2*_sa(r)*∂(r->∂(_Sb,r),r) - p(lb)*r^2*_Sb(r)*∂(r->∂(_sa,r),r))/(r^3*p(lc))
    
    @inline f(i) =  innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTT_pre(lmna, lmnb, lmnc, r,wr, sa,Tb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _sa(i) = sa(la,ma,na,r)
    @inline _Tb(i) = Tb(lb,mb,nb,r)
    @inline _Tc(i) = Tc(lc,mc,nc,r)


    @inline f1(i) = (-p(lc)*(p(la)+p(lb)-p(lc))*(_sa(r)*_Tb(r) + r*∂(_sa,r)*_Tb(r)) +
                        p(la)*(p(la)-p(lb)-p(lc))*(r*∂(_sa,r)*_Tb(r) + r*_sa(r)*∂(_Tb,r)))/(2r^2*p(lc))

    @inline f(i) =  innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end 

#toroidal flow, poloidal B0
function _induction_tST_pre(lmna, lmnb, lmnc, r,wr, ta,Sb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _ta(i) = ta(la,ma,na,r)
    @inline _Sb(i) = Sb(lb,mb,nb,r)
    @inline _Tc(i) = Tc(lc,mc,nc,r)


    @inline f1(i) = (p(lc)*(p(la)+p(lb)-p(lc))*(_ta(r)*_Sb(r) + r*_ta(r)*∂(_Sb,r)) -
                       p(lb)*(p(lb)-p(la)-p(lc))*(r*∂(_ta,r)*_Sb(r) + r*_ta(r)*∂(_Sb,r)))/(2r^2*p(lc))

    @inline f(i) =  innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc

    # aij += f(1.0)*Aabc*lb*lc
    return aij
end 

#toroidal flow, toroidal B0
function _induction_tTT_pre(lmna, lmnb, lmnc, r,wr, ta,Tb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _ta(i) = ta(la,ma,na,r)
    @inline _Tb(i) = Tb(lb,mb,nb,r)
    @inline _Tc(i) = Tc(lc,mc,nc,r)


    @inline f1(i) = _ta(r)*_Tb(r)/r

    @inline f(i) =  innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc

    # aij += f(1.0)*Eabc
    return aij
end 

#create mutating functions for easier sparse matrices assembly
# flist = [:_induction_sSS, :_induction_tSS, :_induction_sTS,
#          :_induction_sST, :_induction_sTT, :_induction_tST, :_induction_tTT]

# for f in flist
#     @eval begin
#         function $(Symbol(string(f)*"!"))(is,js,aijs,i,j, lmna, lmnb, lmnc, r,wr, fa,fb,fc; thresh=sqrt(eps()))
#             aij = $(f)_pre(lmna, lmnb, lmnc, r,wr, fa,fb,fc)
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


# function rhs_induction_bpol_gen(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps())) where T
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = Limace.ChenBasis.lmn_upol(N,m,ns)
#     lmn_t = Limace.ChenBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
#     lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

#     npb = length(lmn_bp)

#     is,js,aijs = Int[],Int[],Complex{T}[]

#     r, wr = rquad(N+lb0+nb0+5)

#     for (i,lmni) in enumerate(lmn_bp)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             # if mj+mb0-mi == 0
#                 _induction_sSS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, s_chen, s_mf, s_mf; thresh)
#             # end
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             _induction_tSS!(is,js,aijs,i,j+np,lmnj,lmnb0,lmni, r, wr, t_chen, s_mf, s_mf; thresh)
#         end
#     end

#     for (i,lmni) in enumerate(lmn_bt)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             # if mj+mb0-mi == 0
#             _induction_sST!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, s_chen, s_mf, t_mf; thresh)
#             # end
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             _induction_tST!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, t_chen, s_mf, t_mf; thresh)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(is,js,aijs,nmatb, nmatu)
# end

# @inline function istriangle(la,lb,lc)
#     (abs(lb-lc)<=la<=lb+lc) && return true
#     return false
# end

# @inline function condition1(la,lb,lc,ma,mb,mc)
#     (mc+mb-ma != 0) && return false
#     isodd(la+lb+lc) && return false
#     !istriangle(la,lb,lc) && return false
#     return true
# end

# @inline function condition1u(la,lb,lc,ma,mb,mc)
#     (mc+ma-mb != 0) && return false
#     isodd(la+lb+lc) && return false
#     !istriangle(la,lb,lc) && return false
#     return true
# end


# @inline function condition2(la,lb,lc,ma,mb,mc)
#     (mc+mb-ma != 0) && return false
#     iseven(la+lb+lc) && return false
#     !istriangle(la,lb,lc) && return false
#     return true
# end

# @inline function condition2u(la,lb,lc,ma,mb,mc)
#     (mc+ma-mb != 0) && return false
#     iseven(la+lb+lc) && return false
#     !istriangle(la,lb,lc) && return false
#     return true
# end


# function ncondition(lb,na,nb,nc)
#     na-nb-lb <= nc <= na+nb+lb
# end

# function _dummy!(is,js,aijs,i,j)
#     push!(is,i)
#     push!(js,j)
#     push!(aijs,1.0)
#     return nothing
# end

function rhs_induction_bpol_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps()), smfb0::Sf = s_mf, conditions=true) where {T,Sf}
    su = s_in
    tu = t_in
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lb0+nb0+5)
    rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
    @inbounds for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            # _dummy!(is,js,aijs,i,j)
            _induction_sSS!(is,js,aijs,i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, s_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            # _dummy!(is,js,aijs,i,j+np)
            _induction_tSS!(is,js,aijs,i,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, s_mf; thresh)
        end
    end

    @inbounds for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sST!(is,js,aijs,i+npb,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            conditions && !ncondition(lb0,ni,nb0,nj) && continue
            conditions && !condition1(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tST!(is,js,aijs,i+npb,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_btor_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps())) where T
    su = s_in
    tu = t_in
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lb0+nb0+5)
    rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, s_mf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTT!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, t_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tTT!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, t_mf, t_mf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_upol_pre(N,m, lmnu0; ns = false, η::T=1.0, thresh = sqrt(eps()), su = s_chen, smf = s_mf, tmf = t_mf, condition = true) where T
    lu0,mu0,nu0 = lmnu0

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    rwrs = [rquad(n+lu0+nu0+1) for n in 1:N]

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sSS!(is,js,aijs,i,j,lmnu0,lmnj,lmni, r, wr, su, smf, smf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTS!(is,js,aijs,i,j+np,lmnu0,lmnj,lmni, r, wr, su, tmf, smf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sST!(is,js,aijs,i+np,j,lmnu0,lmnj,lmni, r, wr, su, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sTT!(is,js,aijs,i+np,j+np,lmnu0,lmnj,lmni, r, wr, su, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end

function rhs_induction_utor_pre(N,m, lmnu0; ns = false, η::T=1.0, thresh = sqrt(eps()), tu = t_chen, smf = s_mf, tmf = t_mf, condition=true) where T
    lu0,mu0,nu0 = lmnu0

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    rwrs = [rquad(n+lu0+nu0+1) for n in 1:N]

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tSS!(is,js,aijs,i,j,lmnu0,lmnj,lmni, r, wr, tu, smf, smf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tST!(is,js,aijs,i+np,j,lmnu0,lmnj,lmni, r, wr, tu, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tTT!(is,js,aijs,i+np,j+np,lmnu0,lmnj,lmni, r, wr, tu, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end


function rhs_induction_btor_cond_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps())) where T
    su = s_in
    tu = t_in
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    tmf = t_in
    smf = s_in
    lmn_bp = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_bt = Limace.InviscidBasis.lmn_utor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _induction_sTS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, su, tmf, smf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _induction_sTT!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, su, tmf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _induction_tTT!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end


function rhs_induction_bpol_dist_pre(N,m, lmnb0; ns = false, η::T=1.0, thresh = sqrt(eps()), smfb0::Sf = s_mf) where {T,Sf}
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
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, s_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j+np)
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, s_mf; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_bt))
        lmni = lmn_bt[i]
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_sST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
            _induction_tST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatb, nmatu)
end

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
