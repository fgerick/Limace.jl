##
## poloidal Lorentz force
##


#poloidal B1, poloidal B2
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_SSs(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Sa = r->s(TA,la,ma,na,r)
    @inline _Sb = r->s(TB,lb,mb,nb,r)
    @inline _sc = r->s(TC,lc,mc,nc,r)

    @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*D(_Sa,la,r)*∂(r->r*_Sb(r),r) + 
                        p(lb)*(p(la)-p(lb)+p(lc))*r*∂(r->D(_Sa,la,r)*_Sb(r),r))/(2r^2*p(lc))

    @inline f = r-> -innert(_sc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end


#poloidal B1, toroidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_STs(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc


    @inline _Sa = r->s(TA,la,ma,na,r)
    @inline _Tb = r->t(TB,lb,mb,nb,r)
    @inline _sc = r->s(TC,lc,mc,nc,r)

    @inline f1 = r -> (p(lc)*r^2*D(_Sa,la,r)*_Tb(r) + 
                        (p(la)+p(lb)+p(lc))*_Sa(r)*_Tb(r) - 
                        (p(la)+p(lb)-p(lc))*(r*_Sa(r)*∂(_Tb,r) + 
                                            r*∂(_Sa,r)*_Tb(r) +
                                            r^2*∂(_Sa,r)*∂(_Tb,r)) - 
                        p(lb)*r^2*∂(r->∂(_Sa,r),r)*_Tb(r) - 
                        p(la)*r^2*∂(r->∂(_Tb,r),r)*_Sa(r)
                        )/(r^3*p(lc))
    
    @inline f = r-> -innert(_sc, f1, lc,r)
    
    aij = ∫dr(f,r,wr)
    return aij
end


#toroidal B1, toroidal B0
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_TTs(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Ta = r->t(TA,la,ma,na,r)
    @inline _Tb = r->t(TB,lb,mb,nb,r)
    @inline _sc = r->s(TC,lc,mc,nc,r)


    @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*∂(r->r*_Ta(r),r)*_Tb(r) + p(la)*(-p(la)+p(lb)+p(lc))*r*∂(r->_Ta(r)*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r->-innert(_sc,f1,lc,r)
    
    aij = ∫dr(f,r,wr)
    return aij
end



##
## toroidal lorentz equation
##


#poloidal B1, poloidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_SSt(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Sa = r->s(TA,la,ma,na,r)
    @inline _Sb = r->s(TB,lb,mb,nb,r)
    @inline _tc = r->t(TC,lc,mc,nc,r)


    @inline f1 = r -> -p(lb)*D(_Sa,la,r)*_Sb(r)/(r*p(lc))
    
    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end


#poloidal B1, toroidal B0
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_STt(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Sa = r->s(TA,la,ma,na,r)
    @inline _Tb = r->t(TB,lb,mb,nb,r)
    @inline _tc = r->t(TC,lc,mc,nc,r)


    @inline f1 = r -> (p(lb)*(p(lb)-p(la)-p(lc))*∂(r->r*_Sa(r),r)*_Tb(r) - p(la)*(p(la)-p(lb)-p(lc))*_Sa(r)*∂(r->r*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end 


#toroidal B1, toroidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_TTt(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Ta = r->t(TA,la,ma,na,r)
    @inline _Tb = r->t(TB,lb,mb,nb,r)
    @inline _tc = r->t(TC,lc,mc,nc,r)


    @inline f1 = r -> p(la)*_Ta(r)*_Tb(r)/(r*p(lc))

    @inline f = r-> innert(_tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end 

#matrix assembly



"""
$(TYPEDSIGNATURES)

Computes the Lorentz term for a poloidal background magnetic field `B0`, a velocity basis `bui` and a magnetic field basis `bbj`.
"""
function lorentz(bui::TI, bbj::TJ, B0::BasisElement{T0,Poloidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    @assert bui.N == bbj.N "Use same resolution for bases!"
    N = bui.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, nrange_p, _lorentz_SSs, A)
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, nrange_p, _lorentz_SSs, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p, nrange_t, _lorentz_STs, E)
        end
    end

    for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, nrange_t, _lorentz_STt, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t, nrange_p, _lorentz_SSt, E)
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t, nrange_p, _lorentz_SSt, E)
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

"""
$(TYPEDSIGNATURES)

Computes the Lorentz term for a toroidal background magnetic field `B0`, a velocity basis `bui` and a magnetic field basis `bbj`.
"""
function lorentz(bui::TI, bbj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    @assert bui.N == bbj.N "Use same resolution for bases!"
    N = bui.N
    rwrs = [rquad(n + l0 + n0 + 5) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p, nrange_t, _lorentz_TTs, A)
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p, nrange_t, _lorentz_TTs, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, nrange_p, _lorentz_STs, E)
        end
    end

    for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t, nrange_p, _lorentz_STt, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, nrange_t, _lorentz_TTt, E)
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, nrange_t, _lorentz_TTt, E)
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end