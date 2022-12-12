@inline p(l) = l*(l+1)

## poloidal induction equation
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

    # @show f(0.5)*(4pi)^(3/2)
    aij = ∫dr(f,r,wr)*Aabc
    return aij
end

function _induction_sSS!(is,js,aijs,i,j, lmna, lmnb, lmnc, r,wr, sa,Sb,Sc; thresh=sqrt(eps()))

    aij = _induction_sSS(lmna, lmnb, lmnc, r,wr, sa,Sb,Sc)
    # if abs(out)>thresh
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    # end
    return nothing
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

function _induction_sTS!(is,js,aijs,i,j, lmna, lmnb, lmnc, r,wr, sa,Tb,Sc)

    aij = _induction_sTS(lmna, lmnb, lmnc, r,wr, sa,Tb,Sc)
    # if abs(out)>thresh
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    # end
    return nothing
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

function _induction_tSS!(is,js,aijs,i,j, lmna, lmnb, lmnc, r,wr, ta,Sb,Sc)

    aij = _induction_tSS(lmna, lmnb, lmnc, r,wr, ta,Sb,Sc)
    # if abs(out)>thresh
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    # end
    return nothing
end

#toroidal flow, toroidal B0 
#always 0
function _induction_tTS!(is,js,aijs,i,j, lmna, lmnb, lmnc)
    return nothing
end


#toroidal induction equation

@inline function innert(t,t2, l, r) 
    return l*(l+1)*t(r)*t2(r)
end


#poloidal flow, poloidal B0
function _induction_sST(lmna, lmnb, lmnc, r,wr, sa,Sb,Tc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = elsasser(la,lb,lc,ma,mb,mc)
    @inline _sa = r->sa(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)
    @inline _Tc = r->Tc(lc,mc,nc,r)


    @inline f1 = r -> ((p(la)+p(lb)+p(lc))*_sa(r)*_Sb(r) - 
                       (p(la)+p(lb)-p(lc))*(r*∂(_sa,r)*_Sb(r) + r*_sa(r)*∂(_Sb,r)+r^2*∂(_sa,r)*∂(_Sb,r)) - p(la)*r^2*_sa(r)*∂(r->∂(_Sb,r),r) - p(lb)*r^2*_Sb(r)*∂(r->∂(_sa,r),r))/(r^3*p(lc))
    
    @inline f = r-> innert(_Tc,f1, lc, r)

    # @show f(0.5)*(4pi)^(3/2)
    aij = ∫dr(f,r,wr)*Aabc
    return aij
end