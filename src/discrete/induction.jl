@inline p(l) = l*(l+1)
## poloidal induction equation
@inline dr(f::T,l,m,n,r) where T = ∂(r->r*f(l,m,n,r),r)

#poloidal flow, poloidal B0
function _induction_sSS(lmna, lmnb, lmnc, r,wr, sa,Sb,Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Aabc = adamgaunt(la,lb,lc,ma,mb,mc)

    @inline f = r -> (-p(la)*(-p(la)+p(lb)+p(lc))*sa(la,ma,na,r)*∂(r->r*Sb(lb,mb,nb,r),r) + 
                      p(lb)*( p(la)-p(lb)+p(lc))*Sb(lb,mb,nb,r)*∂(r->r*sa(la,ma,na,r),r))*Sc(lc,mc,nc,r)/r^2
    
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
function _induction_sTS!(is,js,aijs,i,j, lmna, lmnb, lmnc)

end

#toroidal flow, poloidal B0
function _induction_tSS(lmna, lmnb, lmnc, r,wr, ta,Sb,Sc)
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    Eabc = elsasser(la,lb,lc,ma,mb,mc)

    @inline f = r -> p(lb)*ta(la,ma,na,r)*Sb(lb,mb,nb,r)*Sc(lc,mc,nc,r)/r
    
    aij = ∫dr(f,r,wr)*Eabc
    return aij
end

function _induction_tSS!(is,js,aijs,i,j, lmna, lmnb, lmnc)

end

#toroidal flow, toroidal B0
function _induction_tTS!(is,js,aijs,i,j, lmna, lmnb, lmnc)

end