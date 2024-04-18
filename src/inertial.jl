function _inertial_ss(::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _s = r->s(T,la,ma,na,r)

    @inline f = r-> inners(_s,_s, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end

function _inertial_tt(::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _t = r->t(T,la,ma,na,r)

    @inline f = r-> innert(_t,_t, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end