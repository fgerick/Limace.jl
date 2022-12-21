@inline p(l) = l*(l+1.0)

##
## poloidal induction equation
##

@inline dr(f::T,l,m,n,r) where T = ∂(r->r*f(l,m,n,r),r)

@inline function inners(s,s2, l, r) 
    return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r*s(r),r)*∂(r->r*s2(r),r))/r^2
end


#poloidal flow, poloidal B0
function _induction_sSS(lmna, lmnb, lmnc, r,wr, sa,Sb,Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)

    @inline f1 = r -> (-p(la)*(-p(la)+p(lb)+p(lc))*sa(la,ma,na,r)*∂(r->r*Sb(lb,mb,nb,r),r) + 
                      p(lb)*( p(la)-p(lb)+p(lc))*Sb(lb,mb,nb,r)*∂(r->r*sa(la,ma,na,r),r))/(2r^2*p(lc))
    
    @inline f = r-> inners(r->Sc(lc,mc,nc,r),f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTS(lmna, lmnb, lmnc, r,wr, sa,Tb,Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    @inline f1 = r -> p(la)*sa(la,ma,na,r)*Tb(lb,mb,nb,r)/(r*p(lc))
    @inline f = r->inners(r->Sc(lc,mc,nc,r),f1,lc,r)
    
    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#toroidal flow, poloidal B0
function _induction_tSS(lmna, lmnb, lmnc, r,wr, ta,Sb,Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    @inline f1 = r -> p(lb)*ta(la,ma,na,r)*Sb(lb,mb,nb,r)/(r*p(lc))
    @inline f = r->inners(r->Sc(lc,mc,nc,r),f1,lc,r)
    
    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#toroidal flow, toroidal B0 
#always 0


##
## toroidal induction equation
##

@inline function innert(t,t2, l, r) 
    return l*(l+1)*t(r)*t2(r)
end


#poloidal flow, poloidal B0
function _induction_sST(lmna, lmnb, lmnc, r,wr, sa,Sb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _sa = r->sa(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)
    @inline _Tc = r->Tc(lc,mc,nc,r)


    @inline f1 = r -> ((p(la)+p(lb)+p(lc))*_sa(r)*_Sb(r) - 
                       (p(la)+p(lb)-p(lc))*(r*∂(_sa,r)*_Sb(r) + r*_sa(r)*∂(_Sb,r) + r^2*∂(_sa,r)*∂(_Sb,r)) - p(la)*r^2*_sa(r)*∂(r->∂(_Sb,r),r) - p(lb)*r^2*_Sb(r)*∂(r->∂(_sa,r),r))/(r^3*p(lc))
    
    @inline f = r-> innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTT(lmna, lmnb, lmnc, r,wr, sa,Tb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _sa = r->sa(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _Tc = r->Tc(lc,mc,nc,r)


    @inline f1 = r -> (-p(lc)*(p(la)+p(lb)-p(lc))*(_sa(r)*_Tb(r) + r*∂(_sa,r)*_Tb(r)) +
                        p(la)*(p(la)-p(lb)-p(lc))*(r*∂(_sa,r)*_Tb(r) + r*_sa(r)*∂(_Tb,r)))/(2r^2*p(lc))

    @inline f = r-> innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end 

#toroidal flow, poloidal B0
function _induction_tST(lmna, lmnb, lmnc, r,wr, ta,Sb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _ta = r->ta(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)
    @inline _Tc = r->Tc(lc,mc,nc,r)


    @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*(_ta(r)*_Sb(r) + r*_ta(r)*∂(_Sb,r)) -
                       p(lb)*(p(lb)-p(la)-p(lc))*(r*∂(_ta,r)*_Sb(r) + r*_ta(r)*∂(_Sb,r)))/(2r^2*p(lc))

    @inline f = r-> innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end 

#toroidal flow, toroidal B0
function _induction_tTT(lmna, lmnb, lmnc, r,wr, ta,Tb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _ta = r->ta(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _Tc = r->Tc(lc,mc,nc,r)


    @inline f1 = r -> _ta(r)*_Tb(r)/r

    @inline f = r-> innert(_Tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc
    return aij
end 

#create mutating functions for easier sparse matrices assembly
flist = [:_induction_sSS, :_induction_tSS, :_induction_sTS,
         :_induction_sST, :_induction_sTT, :_induction_tST, :_induction_tTT]

for f in flist
    @eval begin
        function $(Symbol(string(f)*"!"))(is,js,aijs,i,j, lmna, lmnb, lmnc, r,wr, fa,fb,fc; thresh=sqrt(eps()))
            aij = $(f)(lmna, lmnb, lmnc, r,wr, fa,fb,fc)
            if abs(aij) > thresh
                push!(is,i)
                push!(js,j)
                push!(aijs,aij) 
            end
            return nothing
        end
    end
end



#matrix assembly


# function rhs_induction_bpol_gen(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps())) where T
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

@inline function istriangle(la,lb,lc)
    (abs(lb-lc)<=la<=lb+lc) && return true
    return false
end

@inline function condition1(la,lb,lc,ma,mb,mc)
    (mc+mb-ma != 0) && return false
    isodd(la+lb+lc) && return false
    !istriangle(la,lb,lc) && return false
    return true
end

@inline function condition1u(la,lb,lc,ma,mb,mc)
    (mc+ma-mb != 0) && return false
    isodd(la+lb+lc) && return false
    !istriangle(la,lb,lc) && return false
    return true
end


@inline function condition2(la,lb,lc,ma,mb,mc)
    (mc+mb-ma != 0) && return false
    iseven(la+lb+lc) && return false
    !istriangle(la,lb,lc) && return false
    return true
end

@inline function condition2u(la,lb,lc,ma,mb,mc)
    (mc+ma-mb != 0) && return false
    iseven(la+lb+lc) && return false
    !istriangle(la,lb,lc) && return false
    return true
end


function ncondition(lb,na,nb,nc)
    nc <= na+nb+lb
end

function _dummy!(is,js,aijs,i,j)
    push!(is,i)
    push!(js,j)
    push!(aijs,1.0)
    return nothing
end

function rhs_induction_bpol(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps())) where T
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.ChenBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.ChenBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    @inbounds for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _induction_sSS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, s_chen, s_mf, s_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j+np)
            _induction_tSS!(is,js,aijs,i,j+np,lmnj,lmnb0,lmni, r, wr, t_chen, s_mf, s_mf; thresh)
        end
    end

    @inbounds for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _induction_sST!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, s_chen, s_mf, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _induction_tST!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, t_chen, s_mf, t_mf; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_btor(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps())) where T
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.ChenBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.ChenBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _induction_sTS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, s_chen, t_mf, s_mf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _induction_sTT!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, s_chen, t_mf, t_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _induction_tTT!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, t_chen, t_mf, t_mf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatb, nmatu)
end

function rhs_induction_upol(N,m, lmnu0; ns = 0, η::T=1.0, thresh = sqrt(eps()), su = s_chen, smf = s_mf, tmf = t_mf, condition = true) where T
    lu0,mu0,nu0 = lmnu0

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lu0+nu0+5)

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_sSS!(is,js,aijs,i,j,lmnu0,lmnj,lmni, r, wr, su, smf, smf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_sTS!(is,js,aijs,i,j+np,lmnu0,lmnj,lmni, r, wr, su, tmf, smf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_sST!(is,js,aijs,i+np,j,lmnu0,lmnj,lmni, r, wr, su, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_sTT!(is,js,aijs,i+np,j+np,lmnu0,lmnj,lmni, r, wr, su, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end

function rhs_induction_utor(N,m, lmnu0; ns = 0, η::T=1.0, thresh = sqrt(eps()), tu = t_chen, smf = s_mf, tmf = t_mf, condition=true) where T
    lu0,mu0,nu0 = lmnu0

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lu0+nu0+5)

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_tSS!(is,js,aijs,i,j,lmnu0,lmnj,lmni, r, wr, tu, smf, smf; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_tST!(is,js,aijs,i+np,j,lmnu0,lmnj,lmni, r, wr, tu, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0,li,lj,mu0,mi,mj) && continue
            _induction_tTT!(is,js,aijs,i+np,j+np,lmnu0,lmnj,lmni, r, wr, tu, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end