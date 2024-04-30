##
## poloidal induction equation
##

"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(pⱼ×Sₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_sSS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> (-p(la) * (-p(la) + p(lb) + p(lc)) * s(TA, la, ma, na, r) * ∂(r -> r * s(TB, lb, mb, nb, r), r) +
                       p(lb) * (p(la) - p(lb) + p(lc)) * s(TB, lb, mb, nb, r) * ∂(r -> r * s(TA, la, ma, na, r), r)) / (2r^2 * p(lc))

    @inline f = r -> inners(f1, r -> s(TC, lc, mc, nc, r), lc, r)

    aij = ∫dr(f, r, wr) 

    if external
        aij += f1(1.0) * s(TC, lc, mc, nc, 1.0) * p(lc) * lc 
    end
    return aij
end


#poloidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(pⱼ×Tₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_sTS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> p(la) * s(TA, la, ma, na, r) * t(TB, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end


#toroidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(qⱼ×Sₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_tSS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> p(lb) * t(TA, la, ma, na, r) * s(TB, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 

    #add contribution from external ∫dV (1<r<∞), 
    #if toroidal velocity is not 0 at r=11 

    if external
        aij += f1(1.0) * s(TC,lc, mc, nc, 1.0) * lc * p(lc) 
    end

    return aij
end


#toroidal flow, toroidal B0 
#always 0


##
## toroidal induction equation
##

#poloidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(pⱼ×Sₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_sST(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _sa = r -> s(TA, la, ma, na, r)
    @inline _Sb = r -> s(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)

    @inline f1 = r -> ((p(la) + p(lb) + p(lc)) * _sa(r) * _Sb(r) -
                       (p(la) + p(lb) - p(lc)) * (r * ∂(_sa, r) * _Sb(r) + r * _sa(r) * ∂(_Sb, r) + r^2 * ∂(_sa, r) * ∂(_Sb, r)) -
                       p(la) * r^2 * _sa(r) * ∂(r -> ∂(_Sb, r), r) - p(lb) * r^2 * _Sb(r) * ∂(r -> ∂(_sa, r), r)) / (r^3 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end


#poloidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(pⱼ×Tₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_sTT(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _sa = r -> s(TA, la, ma, na, r)
    @inline _Tb = r -> t(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)

    @inline f1 = r -> (-p(lc) * (p(la) + p(lb) - p(lc)) * (_sa(r) * _Tb(r) + r * ∂(_sa, r) * _Tb(r)) +
                       p(la) * (p(la) - p(lb) - p(lc)) * (r * ∂(_sa, r) * _Tb(r) + r * _sa(r) * ∂(_Tb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end

#toroidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(qⱼ×Sₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_tST(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _ta = r -> t(TA, la, ma, na, r)
    @inline _Sb = r -> s(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)

    @inline f1 = r -> (p(lc) * (p(la) + p(lb) - p(lc)) * (_ta(r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r)) -
                       p(lb) * (p(lb) - p(la) - p(lc)) * (r * ∂(_ta, r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end

#toroidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(qⱼ×Tₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_tTT(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _ta = r -> t(TA, la, ma, na, r)
    @inline _Tb = r -> t(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)

    @inline f1 = r -> _ta(r) * _Tb(r) / r
    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end



#matrix assembly

"""
$(TYPEDSIGNATURES)

For a combination of `mi`, `m0` and `mj` only `mj=mi-m0` is nonzero for the Adam-Gaunt variables.
"""
@inline adamgaunt_mjs(mi, m0) = mi - m0

"""
$(TYPEDSIGNATURES)

For a combination of `li`, `l0` and `lj` only the even `lj` that satisfy `|l0-li|≤lj≤l0+li` are nonzero for the Adam-Gaunt variables.
"""
@inline function adamgaunt_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if iseven(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end

"""
$(TYPEDSIGNATURES)

For a combination of `mi`, `m0` and `mj` only `mj=mi-m0` is nonzero for the Elsasser variables.
"""
@inline elsasser_mjs(mi, m0) = mi - m0

"""
$(TYPEDSIGNATURES)

For a combination of `li`, `l0` and `lj` only the odd `lj` that satisfy `|l0-li|≤lj≤l0+li` are nonzero for the Elsasser variables.
"""
@inline function elsasser_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if isodd(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end

"""
$(TYPEDSIGNATURES)

Fallback functions for `_crossterm` term for `U0`. Write specialized function to include e.g. bandedness in `n`.
"""
@inline function _crossterm!(bi::TI, B0::BasisElement{T0,PT,T}, bj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(bi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(T0, TJ, TI, B0.lmn, lmnj, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fallback functions for `_crossterm!` term for `B0`. Write specialized function to include e.g. bandedness in `n`.
"""
@inline function _crossterm!(bi::TI, bj::TJ, B0::BasisElement{T0,PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(bi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end

#∫Sᵢ⋅∇×(pⱼ×S₀) dV
function _induction_sps!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sSS, args...; external)
end

#∫Sᵢ⋅∇×(p₀×Sⱼ) dV
function _induction_sps!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_sSS, args...; external)
end

#∫Sᵢ⋅∇×(pⱼ×T₀) dV
function _induction_spt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sTS, args...)
end

#∫Sᵢ⋅∇×(p₀×Tⱼ) dV
function _induction_spt!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, nrange_p_bc, nrange_t, _induction_sTS, args...)
end

#∫Sᵢ⋅∇×(qⱼ×S₀) dV
function _induction_sqs!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, nrange_p_bc, nrange_t, _induction_tSS, args...; external)
end

#∫Sᵢ⋅∇×(q₀×Sⱼ) dV
function _induction_sqs!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_tSS, args...; external)
end


#∫Tᵢ⋅∇×(pⱼ×S₀) dV
function _induction_tps!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sST, args...)
end

#∫Tᵢ⋅∇×(p₀×Sⱼ) dV
function _induction_tps!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_sST, args...)
end

#∫Tᵢ⋅∇×(pⱼ×T₀) dV
function _induction_tpt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sTT, args...)
end

#∫Tᵢ⋅∇×(p₀×Tⱼ) dV
function _induction_tpt!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_sTT, args...)
end

#∫Tᵢ⋅∇×(qⱼ×S₀) dV
function _induction_tqs!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tST, args...)
end

#∫Tᵢ⋅∇×(q₀×Sⱼ) dV
function _induction_tqs!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_tST, args...)
end

#∫Tᵢ⋅∇×(qⱼ×T₀) dV
function _induction_tqt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tTT, args...)
end

#∫Tᵢ⋅∇×(q₀×Tⱼ) dV
function _induction_tqt!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _crossterm!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_tTT, args...)
end


"""
$(TYPEDSIGNATURES)

Computes the induction term for a poloidal background magnetic field `B0`, a magnetic field basis `bbi` and a velocity basis `buj`.
"""
function induction(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sSS,A; external)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj,nrange_p_bc, nrange_t, _induction_tSS, E; external)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj,nrange_t_bc, nrange_t, _induction_tST, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj,nrange_t_bc, nrange_p, _induction_sST, E)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end

"""
$(TYPEDSIGNATURES)

Computes the induction term for a toroidal background magnetic field `B0`, a magnetic field basis `bbi` and a velocity basis `buj`.
"""
function induction(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sTS, E)
            # _induction_spt!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj,E)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sTT, A)
            # _induction_tpt!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj,A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tTT, E)
            # _induction_tqt!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj,E)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end

"""
$(TYPEDSIGNATURES)

Computes the induction term for a poloidal background velocity `U0`, a magnetic field basis `bbi` and a magnetic field basis `bbj`.
"""
function induction(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _induction_sps!(bbi, U0, bbj, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, A; external)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _induction_spt!(bbi, U0, bbj, is, js, aijs, 0, npbj, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj,E)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _induction_tpt!(bbi, U0, bbj, is, js, aijs, npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj,A)
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(is, js, aijs, nmatbi, nmatbj)
end

"""
$(TYPEDSIGNATURES)

Computes the induction term for a toroidal background velocity `U0`, a magnetic field basis `bbi` and a magnetic field basis `bbj`.
"""
function induction(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _induction_sqs!(bbi, U0, bbj, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, E; external)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _induction_tqs!(bbi, U0, bbj, is, js, aijs, npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj,A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _induction_tqt!(bbi, U0, bbj, is, js, aijs, npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj,E)
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(is, js, aijs, nmatbi, nmatbj)
end

"""
$(TYPEDSIGNATURES)

Threaded version of [induction](@ref)
"""
function induction_threaded(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)

    @sync begin
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_sps!(bbi, buj, B0, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, A; external)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_sqs!(bbi, buj, B0, is[id], js[id], aijs[id], 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, E; external)
            end
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tqs!(bbi, buj, B0, is[id], js[id], aijs[id], npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tps!(bbi, buj, B0, is[id], js[id], aijs[id], npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, E)
            end
        end
    end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
end

"""
$(TYPEDSIGNATURES)

Threaded version of [induction](@ref)
"""
function induction_threaded(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_spt!(bbi, buj, B0, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj,E; external)
            end
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tpt!(bbi, buj, B0, is[id], js[id], aijs[id], npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj,A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tqt!(bbi, buj, B0, is[id], js[id], aijs[id], npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj,E)
            end
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
end

"""
$(TYPEDSIGNATURES)

Threaded version of [induction](@ref)
"""
function induction_threaded(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_sps!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, A; external)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_spt!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, npbj, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj,E; external)
            end
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tpt!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj,A)
            end
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatbi, nmatbj)
end

"""
$(TYPEDSIGNATURES)

Threaded version of [induction](@ref)
"""
function induction_threaded(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_sqs!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, E; external)
            end
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tqs!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj,A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _induction_tqt!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj,E)
            end
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatbi, nmatbj)
end

# function induction_threaded_new(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true, threading=false) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

#     if threading
#         _nt = Threads.nthreads()
#         is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]
#     else
#         is, js, aijs = Int[], Int[], complex(T)[] 
#     end


#     lmn2k_p_bi = lmn2k_p_dict(bbi)
#     lmn2k_t_bi = lmn2k_t_dict(bbi)

#     lmn2k_p_uj = lmn2k_p_dict(buj)
#     lmn2k_t_uj = lmn2k_t_dict(buj)

#     l0, m0, n0 = B0.lmn
#     @assert bbi.N == buj.N
#     N = bbi.N
#     rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]

#     npb = length(lmn2k_p_bi)
#     npu = length(lmn2k_p_uj)

#     @sync begin
#     for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
#         mj = adamgaunt_mjs(mi, m0)
#         for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
#             A = adamgaunt(lj,l0,li, mj, m0, mi)
#             @inline _f = (is,js,aijs) ->  _induction_sps!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, A; external)
#             if threading
#                 Threads.@spawn begin
#                     id = Threads.threadid()
#                    _f(is[id], js[id], aijs[id])
#                 end
#             else
#                 _f(is,js,aijs)
#                 # _induction_sps!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, A; external)
#             end
#         end
#         mj = elsasser_mjs(mi, m0)
#         for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
#             E = elsasser(lj, l0, li, mj, m0, mi)
#             @inline _f = (is,js,aijs) -> _induction_sqs!(bbi, buj, B0, is, js, aijs, 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, E; external)
#             if threading
#                 Threads.@spawn begin
#                     id = Threads.threadid()
#                    _f(is[id], js[id], aijs[id])
#                 end
#             else
#                 _f(is,js,aijs)
#             end
#         end
#     end

#     for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
#         mj = adamgaunt_mjs(mi, m0)
#         for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
#             A = adamgaunt(lj,l0,li, mj, m0, mi)
#             @inline _f = (is,js,aijs) -> _induction_tqs!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, A)
#             if threading
#                 Threads.@spawn begin
#                     id = Threads.threadid()
#                    _f(is[id], js[id], aijs[id])
#                 end
#             else
#                 _f(is,js,aijs)
#             end

#         end
#         mj = elsasser_mjs(mi, m0)
#         for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
#             E = elsasser(lj, l0, li, mj, m0, mi)
#             @inline _f = (is,js,aijs) -> _induction_tps!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, E)
#             if threading
#                 Threads.@spawn begin
#                     id = Threads.threadid()
#                    _f(is[id], js[id], aijs[id])
#                 end
#             else
#                 _f(is,js,aijs)
#             end
#         end
#     end
#     end

#     nmatb = length(bbi)
#     nmatu = length(buj)
    
#     if threading
#         return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
#     else
#         return sparse(is, js, aijs, nmatb, nmatu)
#     end
# end
