##
## poloidal Lorentz force
##

@inline D(f,l,r) = ∂(r->∂(f,r),r) + 2/r * ∂(f,r) - l*(l+1)/r^2 *f(r)
@inline D(f,df,d2f, l,r) = d2f(r) + 2/r * df(r) - l*(l+1)/r^2 *f(r)

@inline function inners2(s,s2, l, r) 
    return l*(l+1)*(s(r)*D(s2,l,r))
end

#poloidal B1, poloidal B2
function _lorentz_SSs(lmna, lmnb, lmnc, r,wr, Sa,Sb,sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _Sa = r->Sa(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)
    @inline _sc = r->sc(lc,mc,nc,r)

    @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*D(_Sa,la,r)*∂(r->r*_Sb(r),r) + 
                        p(lb)*(p(la)-p(lb)+p(lc))*r*∂(r->D(_Sa,la,r)*_Sb(r),r))/(2r^2*p(lc))

    @inline f = r-> -innert(_sc,f1, lc, r)
    # @inline f = r->_sc(r)*f1(r) #inners(_sc,f1,lc,r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end


#poloidal B1, toroidal B0
function _lorentz_STs(lmna, lmnb, lmnc, r,wr, Sa,Tb,sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    @inline _Sa = r->Sa(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _sc = r->sc(lc,mc,nc,r)

    @inline f1 = r -> (p(lc)*r^2*D(_Sa,la,r)*_Tb(r) + 
                        (p(la)+p(lb)+p(lc))*_Sa(r)*_Tb(r) - 
                        (p(la)+p(lb)-p(lc))*(r*_Sa(r)*∂(_Tb,r) + 
                                            r*∂(_Sa,r)*_Tb(r) + 
                                            r^2*∂(_Sa,r)*∂(_Tb,r)) - 
                        p(lb)*r^2*∂(r->∂(_Sa,r),r)*_Tb(r) - 
                        p(la)*r^2*∂(r->∂(_Tb,r),r)*_Sa(r)
                        )/(r^3*p(lc))
    
    @inline f = r-> -innert(_sc, f1, lc,r)
    
    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#toroidal B1, toroidal B0
function _lorentz_TTs(lmna, lmnb, lmnc, r,wr, Ta,Tb,sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _Ta = r->Ta(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _sc = r->sc(lc,mc,nc,r)


    @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*∂(r->r*_Ta(r),r)*_Tb(r) + p(la)*(-p(la)+p(lb)+p(lc))*r*∂(r->_Ta(r)*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r->-innert(_sc,f1,lc,r)
    
    aij = ∫dr(f,r,wr)*Aabc
    return aij
end



##
## toroidal lorentz equation
##


#poloidal B1, poloidal B0
function _lorentz_SSt(lmna, lmnb, lmnc, r,wr, Sa,Sb,tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _Sa = r->Sa(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)
    @inline _tc = r->tc(lc,mc,nc,r)


    @inline f1 = r -> -p(lb)*D(_Sa,la,r)*_Sb(r)/(r*p(lc))
    
    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc
    return aij
end


#poloidal B1, toroidal B0
function _lorentz_STt(lmna, lmnb, lmnc, r,wr, Sa,Tb,tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
    @inline _Sa = r->Sa(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _tc = r->tc(lc,mc,nc,r)


    @inline f1 = r -> (p(lb)*(p(lb)-p(la)-p(lc))*∂(r->r*_Sa(r),r)*_Tb(r) - p(la)*(p(la)-p(lb)-p(lc))*_Sa(r)*∂(r->r*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Aabc
    return aij
end 


#toroidal B1, toroidal B0
function _lorentz_TTt(lmna, lmnb, lmnc, r,wr, Ta,Tb,tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _Ta = r->Ta(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)
    @inline _tc = r->tc(lc,mc,nc,r)


    @inline f1 = r -> p(la)*_Ta(r)*_Tb(r)/(r*p(lc))

    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)*Eabc
    return aij
end 

#create mutating functions for easier sparse matrices assembly
flist = [:_lorentz_SSs, :_lorentz_STs, :_lorentz_TTs,
         :_lorentz_SSt, :_lorentz_STt, :_lorentz_TTt]

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

function rhs_lorentz_bpol(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps()),smfb0::Sf = s_mf) where {T,Sf}
    su = s_in 
    tu = t_in
    smf = s_mf 
    tmf = t_mf 
    
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    @inbounds for (i,lmni) in enumerate(lmn_p)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0+1,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_SSs!(is,js,aijs,i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, smf, smfb0, su; thresh)
            _lorentz_SSs!(is,js,aijs,i,j,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0, smf, su; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_STs!(is,js,aijs,i,j+npb,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0, tmf, su; thresh)
        end
    end

    @inbounds for (i,lmni) in enumerate(lmn_t)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_SSt!(is,js,aijs,i+np,j,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0,smf,tu; thresh)
            _lorentz_SSt!(is,js,aijs,i+np,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, smf,smfb0,tu; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_STt!(is,js,aijs,i+np,j+npb,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0,tmf,tu; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatu, nmatb)
end

# function rhs_lorentz_bpol2(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps()),smfb0::Sf = s_mf) where {T,Sf}
#     su = s_in 
#     tu = t_in
#     smf = s_mf 
#     tmf = t_mf 
    
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
#     lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)
#     lmn_bp_l = Limace.InsulatingMFBasis.lmn_bpol_l(N,m,ns)
#     LP = length(lmn_bp_l)
#     lmn_bt_l = Limace.InsulatingMFBasis.lmn_btor_l(N,m,ns)
#     LT = length(lmn_bt_l)

#     npb = length(lmn_bp)

#     is,js,aijs = Int[],Int[],Complex{T}[]

#     r, wr = rquad(N+lb0+nb0+5)

#     for (i,lmni) in enumerate(lmn_p)
#         li,mi,ni = lmni
#         for lmnbps in view(lmn_bp_l,max(1,li-lb0):min(LP,li+lb0))
#             for (j, lj, mj, nj) in lmnbps 
#                 # lj,mj,nj = lmnj
#                 lmnj = (lj,mj,nj)
#                 # !ncondition(lb0,ni,nb0+1,nj) && continue
#                 # !condition1(li,lb0,lj,mi,mb0,mj) && continue
#                 _dummy!(is,js,aijs,i,j)
#                 # _lorentz_SSs!(is,js,aijs,i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, smf, smfb0, su; thresh)
#                 # _lorentz_SSs!(is,js,aijs,i,j,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0, smf, su; thresh)
#                 #using i,j indices twice for sparse matrix means values are added!
#             end
#         end
#         for lmnbts in view(lmn_bt_l,max(1,li-lb0):min(LT,li+lb0))
#             for (j,lj,mj,nj) in lmnbts
#                 lmnj = (lj,mj,nj)
#                 # !ncondition(lb0,ni,nb0,nj) && continue
#                 # !condition2(li,lb0,lj,mi,mb0,mj) && continue
#                 _dummy!(is,js,aijs,i,j+npb)
#                 # _lorentz_STs!(is,js,aijs,i,j+npb,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0, tmf, su; thresh)
#             end
#         end
#     end

#     for (i,lmni) in enumerate(lmn_t)
#         li,mi,ni = lmni
#         for lmnbps in view(lmn_bp_l,max(1,li-lb0):min(LP,li+lb0))
#             for (j, lj, mj, nj) in lmnbps 
#                 lmnj = (lj,mj,nj)
#                 # !ncondition(lb0,ni,nb0,nj) && continue
#                 # !condition2(li,lb0,lj,mi,mb0,mj) && continue
#                 _dummy!(is,js,aijs,i+np,j)
#                 # _lorentz_SSt!(is,js,aijs,i+np,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, smfb0,smf,tu; thresh)
#                 # _lorentz_SSt!(is,js,aijs,i+np,j,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smf,smfb0,tu; thresh)
#             end
#         end
#         for lmnbts in view(lmn_bt_l,max(1,li-lb0):min(LT,li+lb0))
#             for (j,lj,mj,nj) in lmnbts
#                 lmnj = (lj,mj,nj)
#                 # !ncondition(lb0,ni,nb0,nj) && continue
#                 # !condition1(li,lb0,lj,mi,mb0,mj) && continue
#                 _dummy!(is,js,aijs,i+np,j+npb)
#                 # _lorentz_STt!(is,js,aijs,i+np,j+npb,T.(lmnb0),T.(lmnj),T.(lmni), r, wr, smfb0,tmf,tu; thresh)
#             end
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(is,js,aijs,nmatu, nmatb)
# end


function rhs_lorentz_btor(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps())) where T
    su=s_in 
    tu = t_in
    smf = s_mf 
    tmf = t_mf 
    tmfb0 = t_mf
    
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    for (i,lmni) in enumerate(lmn_p)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_STs!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, smf, tmfb0, su; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_TTs!(is,js,aijs,i,j+npb,lmnj,lmnb0,lmni, r, wr, tmf, tmfb0, su; thresh)
            _lorentz_TTs!(is,js,aijs,i,j+npb,lmnb0,lmnj,lmni, r, wr, tmfb0, tmf, su; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_t)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_STt!(is,js,aijs,i+np,j,lmnj,lmnb0,lmni, r, wr, smf, tmfb0, tu; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_TTt!(is,js,aijs,i+np,j+npb,lmnj,lmnb0,lmni, r, wr, tmf, tmfb0, tu; thresh)
            _lorentz_TTt!(is,js,aijs,i+np,j+npb,lmnb0,lmnj,lmni, r, wr, tmfb0,  tmf, tu; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatu, nmatb)
end


function rhs_lorentz_btor_cond(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps())) where T
    su=s_in 
    tu = t_in
    smf = s_in
    tmf = t_in
    tmfb0 = t_in
    
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_bt = Limace.InviscidBasis.lmn_utor(N,m,ns)

    npb = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    r, wr = rquad(N+lb0+nb0+5)

    for (i,lmni) in enumerate(lmn_p)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_STs!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, smf, tmfb0, su; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_TTs!(is,js,aijs,i,j+npb,lmnj,lmnb0,lmni, r, wr, tmf, tmfb0, su; thresh)
            _lorentz_TTs!(is,js,aijs,i,j+npb,lmnb0,lmnj,lmni, r, wr, tmfb0, tmf, su; thresh)
        end
    end

    for (i,lmni) in enumerate(lmn_t)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_STt!(is,js,aijs,i+np,j,lmnj,lmnb0,lmni, r, wr, smf, tmfb0, tu; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            _lorentz_TTt!(is,js,aijs,i+np,j+npb,lmnj,lmnb0,lmni, r, wr, tmf, tmfb0, tu; thresh)
            _lorentz_TTt!(is,js,aijs,i+np,j+npb,lmnb0,lmnj,lmni, r, wr, tmfb0,  tmf, tu; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(is,js,aijs,nmatu, nmatb)
end


function rhs_lorentz_bpol_dist(N,m, lmnb0; ns = 0, η::T=1.0, thresh = sqrt(eps()),smfb0::Sf = s_mf) where {T,Sf}
    su = s_in 
    tu = t_in
    smf = s_mf 
    tmf = t_mf 
    
    lb0,mb0,nb0 = lmnb0
    lmn_p = Limace.InviscidBasis.lmn_upol(N,m,ns)
    lmn_t = Limace.InviscidBasis.lmn_utor(N,m,ns)

    np = length(lmn_p)

    lmn_bp = Limace.InsulatingMFBasis.lmn_bpol(N,m,ns)
    lmn_bt = Limace.InsulatingMFBasis.lmn_btor(N,m,ns)

    npb = length(lmn_bp)

    nt = nprocs()
    is,js,aijs = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])

    r, wr = rquad(N+lb0+nb0+5)

    @sync @distributed for i in shuffle(eachindex(lmn_p))
        lmni = lmn_p[i]
        it = Threads.threadid()
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0+1,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_SSs!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,lmnj,lmnb0,lmni, r, wr, smf, smfb0, su; thresh)
            _lorentz_SSs!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,lmnb0,lmnj,lmni, r, wr, smfb0, smf, su; thresh)
            #using i,j indices twice for sparse matrix means values are added!
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_STs!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j+npb,lmnb0,lmnj,lmni, r, wr, smfb0, tmf, su; thresh)
        end
    end

    @sync @distributed for i in shuffle(eachindex(lmn_t))
        lmni = lmn_t[i]
        it = Threads.threadid()
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition2(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_SSt!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+np,j,lmnj,lmnb0,lmni, r, wr, smf,smfb0,tu; thresh)
            _lorentz_SSt!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+np,j,lmnb0,lmnj,lmni, r, wr, smfb0,smf,tu; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            !ncondition(lb0,ni,nb0,nj) && continue
            !condition1(li,lb0,lj,mi,mb0,mj) && continue
            # _dummy!(is,js,aijs,i,j)
            _lorentz_STt!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+np,j+npb,lmnb0,lmnj,lmni, r, wr, smfb0,tmf,tu; thresh)
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)
    nmatu = length(lmn_p)+length(lmn_t)

    return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatu, nmatb)
end