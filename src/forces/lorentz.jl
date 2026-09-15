"""
$(TYPEDEF)

- `ubasis::Basis{TU}`: The velocity basis for the Lorentz operator.
- `bbasis::Basis{TB}`: The basis for the Lorentz operator.
- `B0`: The background magnetic field, which can be a single `BasisElement` or a collection of them.
- `factor::T`: : A scalar factor that multiplies the Lorentz operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the Lorentz operator.
- `preassembled::Bool`: A flag indicating whether the Lorentz operator has been preassembled.

## Example usage

```julia
u = Inviscid(10)
b = Insulating(10)
B0 = BasisElement(b, Toroidal, (1,0,1))
f = Limace.Lorentz(u,b,B0)
Limace.assemble!(f)
f.mat = # sparse matrix representation of the Lorentz operator
```

"""
mutable struct Lorentz{TU,TB,T} <: Forcing{2}
    ubasis::Basis{TU}
    bbasis::Basis{TB}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Lorentz(ub::Basis, bb::Basis, B0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(ub), length(bb))
    return Lorentz(ub, bb, B0, ComplexF64(factor), mat, false)
end


function assemble!(f::Lorentz; kwargs...)
    f.mat = sum(lorentz(f.ubasis, f.bbasis, B0; kwargs...) for B0 in f.B0)
    f.preassembled = true
    return f.mat
end

##
## poloidal Lorentz force
##


#poloidal B1, poloidal B2
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_SSs(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline sc = r->s(TC,V,lc,mc,nc,r)
    @inline Sa = r->s(TA,V,la,ma,na,r)
    @inline Sb = r->s(TB,V,lb,mb,nb,r)

    @inline function f1(r)
        _Sa, _dSa, _d2Sa, _d3Sa = derivatives0123(Sa,r)
        _Sb, _dSb = derivatives01(Sb,r)

    
        _Da = D(_Sa,_dSa, _d2Sa, la,r)
        _dDa = dD(_Sa,_dSa, _d2Sa, _d3Sa, la,r)

        return  (p(lc)*(p(la)+p(lb)-p(lc))*_Da*(r*_dSb+_Sb) + 
                            p(lb)*(p(la)-p(lb)+p(lc))*r*(_dDa*_Sb+_Da*_dSb))/(2r^2*p(lc))
    end

    # @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*D(_Sa,la,r)*∂(r->r*_Sb(r),r) + 
    #                     p(lb)*(p(la)-p(lb)+p(lc))*r*∂(r->D(_Sa,la,r)*_Sb(r),r))/(2r^2*p(lc))

    @inline f = r-> -innert(sc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end


#poloidal B1, toroidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_STs(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc


    @inline Sa = r->s(TA,V,la,ma,na,r)
    @inline Tb = r->t(TB,V,lb,mb,nb,r)
    @inline sc = r->s(TC,V,lc,mc,nc,r)

    @inline function f1(r)

        _Sa, _dSa, _d2Sa = derivatives012(Sa,r)
        _Tb, _dTb, _d2Tb = derivatives012(Tb,r)
        _Da = D(_Sa, _dSa, _d2Sa, la, r)
        return (p(lc)*r^2*_Da*_Tb + 
                        (p(la)+p(lb)+p(lc))*_Sa*_Tb - 
                        (p(la)+p(lb)-p(lc))*(r*_Sa*_dTb + 
                                            r*_dSa*_Tb +
                                            r^2*_dSa*_dTb) - 
                        p(lb)*r^2*_d2Sa*_Tb - 
                        p(la)*r^2*_d2Tb*_Sa
                        )/(r^3*p(lc))
    end
    # @inline f1 = r -> (p(lc)*r^2*D(_Sa,la,r)*_Tb(r) + 
    #                     (p(la)+p(lb)+p(lc))*_Sa(r)*_Tb(r) - 
    #                     (p(la)+p(lb)-p(lc))*(r*_Sa(r)*∂(_Tb,r) + 
    #                                         r*∂(_Sa,r)*_Tb(r) +
    #                                         r^2*∂(_Sa,r)*∂(_Tb,r)) - 
    #                     p(lb)*r^2*∂(r->∂(_Sa,r),r)*_Tb(r) - 
    #                     p(la)*r^2*∂(r->∂(_Tb,r),r)*_Sa(r)
    #                     )/(r^3*p(lc))
    
    @inline f = r-> -innert(sc, f1, lc,r)
    
    aij = ∫dr(f,r,wr)
    return aij
end


#toroidal B1, toroidal B0
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_TTs(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline Ta = r->t(TA,V,la,ma,na,r)
    @inline Tb = r->t(TB,V,lb,mb,nb,r)
    @inline sc = r->s(TC,V,lc,mc,nc,r)

    @inline function f1(r)
        _Ta, _dTa = derivatives01(Ta,r)
        _Tb, _dTb = derivatives01(Tb,r)
        return (p(lc)*(p(la)+p(lb)-p(lc))*(r*_dTa+_Ta)*_Tb + p(la)*(-p(la)+p(lb)+p(lc))*r*(_dTa*_Tb+_Ta*_dTb))/(2r^2*p(lc))
    end


    # @inline f1 = r -> (p(lc)*(p(la)+p(lb)-p(lc))*∂(r->r*_Ta(r),r)*_Tb(r) + p(la)*(-p(la)+p(lb)+p(lc))*r*∂(r->_Ta(r)*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r->-innert(sc,f1,lc,r)
    
    aij = ∫dr(f,r,wr)
    return aij
end



##
## toroidal lorentz equation
##


#poloidal B1, poloidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_SSt(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline Sa = r->s(TA,V,la,ma,na,r)
    @inline Sb = r->s(TB,V,lb,mb,nb,r)
    @inline tc = r->t(TC,V,lc,mc,nc,r)


    @inline f1 = r -> -p(lb)*D(Sa,la,r)*Sb(r)/(r*p(lc))
    
    @inline f = r-> innert(tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end


#poloidal B1, toroidal B0
# Aabc = adamgaunt(la,lb,lc,ma,mb,mc)
function _lorentz_STt(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline Sa = r->s(TA,V,la,ma,na,r)
    @inline Tb = r->t(TB,V,lb,mb,nb,r)
    @inline tc = r->t(TC,V,lc,mc,nc,r)

    @inline function f1(r)
        _Sa, _dSa = derivatives01(Sa,r)
        _Tb, _dTb = derivatives01(Tb,r)
        return (p(lb)*(p(lb)-p(la)-p(lc))*(r*_dSa+_Sa)*_Tb - p(la)*(p(la)-p(lb)-p(lc))*_Sa*(r*_dTb+_Tb))/(2r^2*p(lc))
    end

    # @inline f1 = r -> (p(lb)*(p(lb)-p(la)-p(lc))*∂(r->r*_Sa(r),r)*_Tb(r) - p(la)*(p(la)-p(lb)-p(lc))*_Sa(r)*∂(r->r*_Tb(r),r))/(2r^2*p(lc))

    @inline f = r-> innert(tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end 


#toroidal B1, toroidal B0
# Eabc = elsasser(la,lb,lc,ma,mb,mc)
function _lorentz_TTt(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline Ta = r->t(TA,V,la,ma,na,r)
    @inline Tb = r->t(TB,V,lb,mb,nb,r)
    @inline tc = r->t(TC,V,lc,mc,nc,r)


    @inline f1 = r -> p(la)*Ta(r)*Tb(r)/(r*p(lc))

    @inline f = r-> innert(tc,f1, lc, r)

    aij = ∫dr(f,r,wr)
    return aij
end 

#matrix assembly



"""
$(TYPEDSIGNATURES)

Computes the Lorentz term for a poloidal background magnetic field `B0`, a velocity basis `bui` and a magnetic field basis `bbj`.
"""
function _lorentz(::Val{false}, bui::TI, bbj::TJ, B0::BasisElement{T0,Poloidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_SSs, A)
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_SSs, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_STs, E)
        end
    end

    for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_STt, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_SSt, E)
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_SSt, E)
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function _lorentz_new(::Val{false}, bui::TI, bbj::TJ, B0::BasisElement{T0,Poloidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)


    for li in 1:lpmax(bui)
        for ni in nrange_p_bc(bui, li)
            for lj in adamgaunt_ljs(li, l0, 0, lpmax(bbj))
               _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSs)  
               _crossterm_m_adamgaunt!(bui,B0,bbj, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSs)  
            end
            for lj in elsasser_ljs(li, l0, 0, ltmax(bbj))
               _crossterm_m_elsasser!(bui,B0,bbj, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_STs)  
            end
        end
    end

    for li in 1:ltmax(bui)
        for ni in nrange_t_bc(bui, li)
            for lj in  elsasser_ljs(li, l0, 0, lpmax(bbj))
               _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSt)  
               _crossterm_m_elsasser!(bui,B0,bbj, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSt)  
            end
            for lj in  adamgaunt_ljs(li, l0, 0, ltmax(bbj))
               _crossterm_m_adamgaunt!(bui, B0,bbj, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_STt)
            end
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
function _lorentz(::Val{false}, bui::TI, bbj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_TTs, A)
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_TTs, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_STs, E)
        end
    end

    for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_STt, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bui, bbj, B0, is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_TTt, E)
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bui, B0, bbj, is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_TTt, E)
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function _lorentz_new(::Val{false}, bui::TI, bbj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)


    for li in 1:lpmax(bui)
        for ni in nrange_p_bc(bui, li)
            for lj in elsasser_ljs(li, l0, 0, lpmax(bbj))
               _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_STs)  
            end
            for lj in adamgaunt_ljs(li, l0, 0, ltmax(bbj))
               _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTs)  
               _crossterm_m_adamgaunt!(bui,B0,bbj, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTs)  
            end
        end
    end

    for li in 1:ltmax(bui)
        for ni in nrange_t_bc(bui, li)
            for lj in  adamgaunt_ljs(li, l0, 0, lpmax(bbj))
               _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_STt)  
            end
            for lj in  elsasser_ljs(li, l0, 0, ltmax(bbj))
               _crossterm_m_elsasser!(bui, B0,bbj, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTt)
               _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTt)
            end
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

"""
$(TYPEDSIGNATURES)

Computes the Lorentz term for a poloidal background magnetic field `B0`, a velocity basis `bui` and a magnetic field basis `bbj`.
"""
function _lorentz(::Val{true}, bui::TI, bbj::TJ, B0::BasisElement{T0,Poloidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    @sync for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _crossterm!(bui, bbj, B0 , is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_SSs, A)
            end
            Threads.@spawn begin
				A = adamgaunt(l0,lj,li, m0, mj, mi)
                _crossterm!(bui, B0, bbj , is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_SSs, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            Threads.@spawn begin
                E = elsasser(l0, lj, li, m0, mj, mi)
                _crossterm!(bui, B0, bbj , is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_STs, E)
            end
        end
    end

    @sync for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            Threads.@spawn begin
                A = adamgaunt(l0,lj,li, m0, mj, mi)
                _crossterm!(bui, B0, bbj , is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_STt, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bui, bbj, B0 , is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_SSt, E)
            end
            Threads.@spawn begin
                E = elsasser(l0, lj, li, m0, mj, mi)
                _crossterm!(bui, B0, bbj , is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_SSt, E)
            end
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function _lorentz_new(::Val{true}, bui::TI, bbj::TJ, B0::BasisElement{T0,Poloidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)


    @sync begin
        for li in 1:lpmax(bui)
            for ni in nrange_p_bc(bui, li)
                for lj in adamgaunt_ljs(li, l0, 0, lpmax(bbj))
                    Threads.@spawn _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSs)  
                    Threads.@spawn _crossterm_m_adamgaunt!(bui,B0,bbj, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSs)  
                end
                for lj in elsasser_ljs(li, l0, 0, ltmax(bbj))
                    Threads.@spawn _crossterm_m_elsasser!(bui,B0,bbj, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_STs)  
                end
            end
        end

        for li in 1:ltmax(bui)
            for ni in nrange_t_bc(bui, li)
                for lj in  elsasser_ljs(li, l0, 0, lpmax(bbj))
                    Threads.@spawn _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSt)  
                    Threads.@spawn _crossterm_m_elsasser!(bui,B0,bbj, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_SSt)  
                end
                for lj in  adamgaunt_ljs(li, l0, 0, ltmax(bbj))
                    Threads.@spawn _crossterm_m_adamgaunt!(bui, B0,bbj, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_STt)
                end
            end
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
function _lorentz(::Val{true}, bui::TI, bbj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    @sync for li in 1:lpmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _crossterm!(bui, bbj, B0, is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_TTs, A)
            end
            Threads.@spawn begin
                A = adamgaunt(l0,lj,li, m0, mj, mi)
                _crossterm!(bui, B0, bbj, is, js, aijs, lck, 0, npb, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_p_bc, nrange_t, _lorentz_TTs, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bui, bbj, B0, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p_bc, nrange_p, _lorentz_STs, E)
            end
        end
    end

    @sync for li in 1:ltmax(bui), mi in intersect(bui.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _crossterm!(bui, bbj, B0, is, js, aijs, lck, npu, 0, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_t_bc, nrange_p, _lorentz_STt, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bui, bbj, B0, is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_TTt, E)
            end
            Threads.@spawn begin
                E = elsasser(l0, lj, li, m0, mj, mi)
                _crossterm!(bui, B0, bbj, is, js, aijs, lck, npu, npb, li, mi, lj, mj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t_bc, nrange_t, _lorentz_TTt, E)
            end
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

function _lorentz_new(::Val{true}, bui::TI, bbj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(bui)
    lmn2k_t_ui = lmn2k_t_dict(bui)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = B0.lmn
    # @assert bui.N == bbj.N "Use same resolution for bases!"
    N = max(bui.N,bbj.N)
    rwrs = [rquad(n + l0 + n0 + 1, bui.V) for n in 1:N]

    npu = length(lmn2k_p_ui)
    npb = length(lmn2k_p_bj)

    @sync begin
        for li in 1:lpmax(bui)
            for ni in nrange_p_bc(bui, li)
                for lj in elsasser_ljs(li, l0, 0, lpmax(bbj))
                    Threads.@spawn _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_STs)  
                end
                for lj in adamgaunt_ljs(li, l0, 0, ltmax(bbj))
                    Threads.@spawn _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTs)  
                    Threads.@spawn _crossterm_m_adamgaunt!(bui,B0,bbj, is, js, aijs, lck, 0,npb, li,ni,lj, rwrs, lmn2k_p_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTs)  
                end
            end
        end

        for li in 1:ltmax(bui)
            for ni in nrange_t_bc(bui, li)
                for lj in  adamgaunt_ljs(li, l0, 0, lpmax(bbj))
                    Threads.@spawn _crossterm_m_adamgaunt!(bui,bbj,B0, is, js, aijs, lck, npu,0, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_p_bj, nrange_p, lpmax, _lorentz_STt)  
                end
                for lj in  elsasser_ljs(li, l0, 0, ltmax(bbj))
                    Threads.@spawn _crossterm_m_elsasser!(bui, B0,bbj, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTt)
                    Threads.@spawn _crossterm_m_elsasser!(bui,bbj,B0, is, js, aijs, lck, npu, npb, li,ni,lj, rwrs, lmn2k_t_ui, lmn2k_t_bj, nrange_t, ltmax, _lorentz_TTt)
                end
            end
        end
    end

    nmatu = length(bui)
    nmatb = length(bbj)

    return sparse(is, js, aijs, nmatu, nmatb)
end

"""
$(TYPEDSIGNATURES)

Computes the Lorentz term for a poloidal/toroidal background magnetic field `B0`, a velocity basis `bui` and a magnetic field basis `bbj`.
"""
function lorentz(bui::TI, bbj::TJ, B0::BasisElement{T0,TH,T}; threads=false, external=false) where {TI<:Basis,TJ<:Basis,T0<:Basis,TH<:Helmholtz,T}
    if length(bui.m) == length(bbj.m) == 1
        return _lorentz(Val(threads), bui, bbj, B0)
    else
        return _lorentz_new(Val(threads), bui, bbj, B0)
    end
end

