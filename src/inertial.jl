
@inline p(l) = l*(l+1.0)

@inline function innert(t,t2, l, r) 
    return l*(l+1)*t(r)*t2(r)
end

@inline function inners(s,s2, l, r) 
    return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r*s(r),r)*∂(r->r*s2(r),r))/r^2
end

function _inertial_ss(::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(T,la,ma,na,r)
    @inline _sb = r->s(T,lb,mb,nb,r)

    @inline f = r-> inners(_sa,_sb, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end

function _inertial_tt(::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(T,la,ma,na,r)
    @inline _tb = r->t(T,lb,mb,nb,r)

    @inline f = r-> innert(_ta,_tb, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end