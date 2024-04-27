##
## poloidal induction equation
##

#poloidal flow, poloidal B0
function _induction_sSS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Aabc = adamgaunt(la, lb, lc, ma, mb, mc)

    @inline f1 = r -> (-p(la) * (-p(la) + p(lb) + p(lc)) * s(TA, la, ma, na, r) * ∂(r -> r * s(TB, lb, mb, nb, r), r) +
                       p(lb) * (p(la) - p(lb) + p(lc)) * s(TB, lb, mb, nb, r) * ∂(r -> r * s(TA, la, ma, na, r), r)) / (2r^2 * p(lc))

    @inline f = r -> inners(f1, r -> s(TC, lc, mc, nc, r), lc, r)

    aij = ∫dr(f, r, wr) 
    # * Aabc

    if external
        aij += f1(1.0) * s(TC, lc, mc, nc, 1.0) * p(lc) * lc 
        # * Aabc
    end
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Eabc = elsasser(la, lb, lc, ma, mb, mc)

    @inline f1 = r -> p(la) * s(TA, la, ma, na, r) * t(TB, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Eabc
    return aij
end


#toroidal flow, poloidal B0
function _induction_tSS(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Eabc = elsasser(la, lb, lc, ma, mb, mc)

    @inline f1 = r -> p(lb) * t(TA, la, ma, na, r) * s(TB, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Eabc

    #add contribution from external ∫dV (1<r<∞), 
    #if toroidal velocity is not 0 at r=11 

    if external
        # aij += lc*p(lb)*ta(la,ma,na,1.0)*Sb(lb,mb,nb,1.0)*Sc(lc,mc,nc,1.0)*Eabc
        aij += f1(1.0) * s(TC,lc, mc, nc, 1.0) * lc * p(lc) 
        # * Eabc
    end

    return aij
end


#toroidal flow, toroidal B0 
#always 0


##
## toroidal induction equation
##

#poloidal flow, poloidal B0
function _induction_sST(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Eabc = elsasser(la, lb, lc, ma, mb, mc)
    @inline _sa = r -> s(TA, la, ma, na, r)
    @inline _Sb = r -> s(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)


    @inline f1 = r -> ((p(la) + p(lb) + p(lc)) * _sa(r) * _Sb(r) -
                       (p(la) + p(lb) - p(lc)) * (r * ∂(_sa, r) * _Sb(r) + r * _sa(r) * ∂(_Sb, r) + r^2 * ∂(_sa, r) * ∂(_Sb, r)) -
                       p(la) * r^2 * _sa(r) * ∂(r -> ∂(_Sb, r), r) - p(lb) * r^2 * _Sb(r) * ∂(r -> ∂(_sa, r), r)) / (r^3 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Eabc
    return aij
end


#poloidal flow, toroidal B0
function _induction_sTT(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Aabc = adamgaunt(la, lb, lc, ma, mb, mc)
    @inline _sa = r -> s(TA, la, ma, na, r)
    @inline _Tb = r -> t(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)


    @inline f1 = r -> (-p(lc) * (p(la) + p(lb) - p(lc)) * (_sa(r) * _Tb(r) + r * ∂(_sa, r) * _Tb(r)) +
                       p(la) * (p(la) - p(lb) - p(lc)) * (r * ∂(_sa, r) * _Tb(r) + r * _sa(r) * ∂(_Tb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Aabc
    return aij
end

#toroidal flow, poloidal B0
function _induction_tST(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Aabc = adamgaunt(la, lb, lc, ma, mb, mc)
    @inline _ta = r -> t(TA, la, ma, na, r)
    @inline _Sb = r -> s(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)


    @inline f1 = r -> (p(lc) * (p(la) + p(lb) - p(lc)) * (_ta(r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r)) -
                       p(lb) * (p(lb) - p(la) - p(lc)) * (r * ∂(_ta, r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Aabc

    # aij += f(1.0)*Aabc*lb*lc
    return aij
end

#toroidal flow, toroidal B0
function _induction_tTT(::Type{TA}, ::Type{TB}, ::Type{TC}, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    # Eabc = elsasser(la, lb, lc, ma, mb, mc)
    @inline _ta = r -> t(TA, la, ma, na, r)
    @inline _Tb = r -> t(TB, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, lc, mc, nc, r)


    @inline f1 = r -> _ta(r) * _Tb(r) / r

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    # * Eabc

    # aij += f(1.0)*Eabc
    return aij
end








#matrix assembly

@inline adamgaunt_mjs(mi, m0) = mi - m0

@inline function adamgaunt_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if iseven(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end

@inline elsasser_mjs(mi, m0) = mi - m0

@inline function elsasser_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if isodd(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end


@inline function __induction!(bbi::TI, U0::BasisElement{T0,PT,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    for ni in nrangefi(bbi, li)
        for nj in nrangefj(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr; kwargs...)
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij)
        end
    end
    return nothing
end

@inline function __induction!(bbi::TI, bbj::TJ, B0::BasisElement{T0,PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_uj, nrangefi, nrangefj, indf; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    for ni in nrangefi(bbi, li)
        for nj in nrangefj(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr; kwargs...)
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_uj[lmnj] + j0, aij)
        end
    end
    return nothing
end

@inline function __induction!(bbi::TI, bbj::TJ, B0::BasisElement{T0,PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_uj, nrangefi, nrangefj, indf::TF, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T, TF}
    for ni in nrangefi(bbi, li)
        for nj in nrangefj(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_uj[lmnj] + j0, aij)
        end
    end
    return nothing
end

#∫Sᵢ⋅∇×(pⱼ×S₀) dV
function _induction_sps!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p, nrange_p, _induction_sSS, args...; external)
end

#∫Sᵢ⋅∇×(p₀×Sⱼ) dV
function _induction_sps!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p, nrange_p, _induction_sSS, args...; external)
end

#∫Sᵢ⋅∇×(pⱼ×T₀) dV
function _induction_spt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p, nrange_p, _induction_sTS, args...)
end

#∫Sᵢ⋅∇×(p₀×Tⱼ) dV
function _induction_spt!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, nrange_p, nrange_t, _induction_sTS, args...)
end

#∫Sᵢ⋅∇×(qⱼ×S₀) dV
function _induction_sqs!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, nrange_p, nrange_t, _induction_tSS, args...; external)
end

#∫Sᵢ⋅∇×(q₀×Sⱼ) dV
function _induction_sqs!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, args...; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p, nrange_p, _induction_tSS, args...; external)
end


#∫Tᵢ⋅∇×(pⱼ×S₀) dV
function _induction_tps!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t, nrange_p, _induction_sST, args...)
end

#∫Tᵢ⋅∇×(p₀×Sⱼ) dV
function _induction_tps!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t, nrange_p, _induction_sST, args...)
end

#∫Tᵢ⋅∇×(pⱼ×T₀) dV
function _induction_tpt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t, nrange_p, _induction_sTT, args...)
end

#∫Tᵢ⋅∇×(p₀×Tⱼ) dV
function _induction_tpt!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t, nrange_t, _induction_sTT, args...)
end

#∫Tᵢ⋅∇×(qⱼ×S₀) dV
function _induction_tqs!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t, nrange_t, _induction_tST, args...)
end

#∫Tᵢ⋅∇×(q₀×Sⱼ) dV
function _induction_tpt!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t, nrange_p, _induction_tST, args...)
end

#∫Tᵢ⋅∇×(qⱼ×T₀) dV
function _induction_tqt!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t, nrange_t, _induction_tTT, args...)
end

#∫Tᵢ⋅∇×(q₀×Tⱼ) dV
function _induction_tqt!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, args...) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return __induction!(bbi, U0, bbj, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t, nrange_t, _induction_tTT, args...)
end



# function _induction_adamgaunt(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
#     is, js, aijs = Int[], Int[], Complex{T}[]

#     lmn2k_p_bi = lmn2k_p_dict(bbi)
#     lmn2k_p_uj = lmn2k_p_dict(buj)

#     l0, m0, n0 = B0.lmn
#     @assert bbi.N == buj.N
#     N = bbi.N
#     rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]
#     # r, wr = rquad(N + l0+n0+1)

#     nbi = length(lmn2k_p_bi)
#     nuj = length(lmn2k_p_uj)

#     for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
#         mj = adamgaunt_mjs(mi, m0)
#         for lj in adamgaunt_ljs(li, l0, mj)
#             _induction_poloidal_poloidal_poloidal!(bbi, buj, B0, is, js, aijs, i0, j0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj; external)
#         end
#     end

#     return sparse(is, js, aijs, nbi, nuj)
# end

function induction(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], Complex{T}[]

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
    k = 1
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _induction_sps!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj,A; external)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _induction_sqs!(bbi, buj, B0, is, js, aijs, 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj,E; external)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _induction_tqs!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj,A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _induction_tps!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj,E)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end

function induction_threaded(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [Complex{T}[] for _ in 1:_nt]

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

function induction_threaded2(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [Complex{T}[] for _ in 1:_nt]

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

    # @sync begin
    Threads.@threads :static for li in 1:lpmax(bbi)
        id = Threads.threadid()
        for mi in intersect(bbi.m, -li:li)
            mj = adamgaunt_mjs(mi, m0)
            for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _induction_sps!(bbi, buj, B0, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, A; external)
            end
            mj = elsasser_mjs(mi, m0)
            for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
                E = elsasser(lj, l0, li, mj, m0, mi)
                _induction_sqs!(bbi, buj, B0, is[id], js[id], aijs[id], 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, E; external)
            end
        end
    end

    Threads.@threads :static for li in 1:ltmax(bbi)
        id = Threads.threadid()
        for mi in intersect(bbi.m, -li:li)
            mj = adamgaunt_mjs(mi, m0)
            for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _induction_tqs!(bbi, buj, B0, is[id], js[id], aijs[id], npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, A)
            end
            mj = elsasser_mjs(mi, m0)
            for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
                E = elsasser(lj, l0, li, mj, m0, mi)
                _induction_tps!(bbi, buj, B0, is[id], js[id], aijs[id], npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, E)
            end
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
end



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
    na-nb-lb <= nc <= na+nb+lb
end


#∫Sᵢ⋅∇×(pⱼ×S₀) dV
function _induction_poloidal_poloidal_poloidal_classic(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true, conditions=true, nconditions=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmnp_bi = lmn_p(bbi)
    lmnp_uj = lmn_p(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    r, wr = rquad(bbi.N + l0 + n0 + 1)

    nbi = length(lmnp_bi)
    nuj = length(lmnp_uj)

    @inbounds for (i, lmni) in enumerate(lmnp_bi)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmnp_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(l0, ni, n0 + 1, nj) && continue
            conditions && !condition1(li, l0, lj, mi, m0, mj) && continue
            aij = _induction_sSS(TJ, T0, TI, (lj, mj, nj), (l0, m0, n0), (li, mi, ni), r, wr; external)
            appendit!(is, js, aijs, i, j, aij)
        end
    end

    return sparse(is, js, aijs, nbi, nuj)
end

function _induction_poloidal_poloidal_poloidal_classic(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true, conditions=true, nconditions=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmnp_bi = lmn_p(bbi)
    lmnp_bj = lmn_p(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N
    r, wr = rquad(bbi.N + l0 + n0 + 1)

    nbi = length(lmnp_bi)
    nbj = length(lmnp_bj)

    @inbounds for (i, lmni) in enumerate(lmnp_bi)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmnp_bj)
            lj, mj, nj = lmnj
            # conditions && nconditions && !ncondition(l0,ni,n0+1,nj) && continue
            conditions && !condition1u(l0, li, lj, m0, mi, mj) && continue
            aij = _induction_sSS(T0, TJ, TI, (l0, m0, n0), (lj, mj, nj), (li, mi, ni), r, wr; external)
            appendit!(is, js, aijs, i, j, aij)
        end
    end

    return sparse(is, js, aijs, nbi, nbj)
end

#∫Sᵢ⋅∇×(pⱼ×S₀) dV
function _induction_poloidal_poloidal_toroidal_classic(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, conditions=true, nconditions=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmnp_bi = lmn_p(bbi)
    lmnp_uj = lmn_p(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    r, wr = rquad(bbi.N + l0 + n0 + 1)

    nbi = length(lmnp_bi)
    nuj = length(lmnp_uj)

    @inbounds for (i, lmni) in enumerate(lmnp_bi)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmnp_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(l0, ni, n0 + 1, nj) && continue
            conditions && !condition2(li, l0, lj, mi, m0, mj) && continue
            aij = _induction_sTS(TJ, T0, TI, (lj, mj, nj), (l0, m0, n0), (li, mi, ni), r, wr)
            appendit!(is, js, aijs, i, j, aij)
        end
    end

    return sparse(is, js, aijs, nbi, nuj)
end






function induction_classic(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true, conditions=true, nconditions=true, thresh=sqrt(eps())) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], Complex{T}[]

    lmnp_bi = lmn_p(bbi)
    lmnp_uj = lmn_p(buj)

    lmnt_bi = lmn_t(bbi)
    lmnt_uj = lmn_t(buj)

    lb0, mb0, nb0 = lmnb0 = B0.lmn
    @assert bbi.N == buj.N
    N = bbi.N
    rwrs = [rquad(n + lb0 + nb0 + 1) for n in 1:N]

    npb = length(lmnp_bi)
    npu = length(lmnp_uj)

    @inbounds for (i, lmni) in enumerate(lmnp_bi)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmnp_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _dummy!(is,js,aijs,i,j)
            # _induction_sSS!(is,js,aijs,i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, sb; thresh, external)
            # aij = _induction_sSS(T.(lmnj), T.(lmnb0), T.(lmni), r, wr, su, smfb0, sb; external)
            A = adamgaunt(lj,lb0,li, mj, mb0, mi)
            aij = _induction_sSS(TJ,T0,TI,(lj,mj,nj), lmnb0,  (li,mi,ni), r,wr; external)*A
            appendit!(is, js, aijs, i, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmnt_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _dummy!(is,js,aijs,i,j+np)
            # _induction_tSS!(is,js,aijs,i,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, s_mf; thresh, external)
            # aij = _induction_tSS(T.(lmnj), T.(lmnb0), T.(lmni), r, wr, tu, smfb0, sb; external)
            E = elsasser(lj,lb0,li, mj, mb0, mi)
            aij = _induction_tSS(TJ,T0,TI, (lj,mj,nj), lmnb0, (li,mi,ni), r,wr; external)*E
            appendit!(is, js, aijs, i, j + npu, aij; thresh)
        end
    end

    @inbounds for (i, lmni) in enumerate(lmnt_bi)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmnp_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _induction_sST!(is, js, aijs, i + npb, j, T.(lmnj), T.(lmnb0), T.(lmni), r, wr, su, smfb0, tb; thresh)
            E = elsasser(lj,lb0,li, mj, mb0, mi)
            aij = _induction_sST(TJ,T0,TI, (lj,mj,nj), lmnb0, (li,mi,ni), r,wr)*E
            appendit!(is, js, aijs, i + npb, j, aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmnt_uj)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _induction_tST!(is, js, aijs, i + npb, j + npu, T.(lmnj), T.(lmnb0), T.(lmni), r, wr, tu, smfb0, tb; thresh)
            A = adamgaunt(lj,lb0,li, mj, mb0, mi)
            aij = _induction_tST(TJ,T0,TI, (lj,mj,nj), lmnb0, (li,mi,ni), r,wr)*A
            appendit!(is, js, aijs, i + npb, j + npu, aij; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end



function rhs_induction_bpol(N, m, lmnb0; ns=false, thresh::T=sqrt(eps()), smfb0::Sf=s_mf, conditions=true, nconditions=true, external=true,
    su=s_in,
    tu=t_in,
    sb=s_mf,
    tb=t_mf,
    lmn_p=InviscidBasis.lmn_upol(N, m, ns),
    lmn_t=InviscidBasis.lmn_utor(N, m, ns),
    lmn_bp=InsulatingMFBasis.lmn_bpol(N, m, ns),
    lmn_bt=InsulatingMFBasis.lmn_btor(N, m, ns)) where {T,Sf}

    lb0, mb0, nb0 = lmnb0
    np = length(lmn_p)
    npb = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]

    # r, wr = rquad(N+lb0+nb0+5)
    rwrs = [rquad(n + lb0 + nb0 + 1) for n in 1:N]
    @inbounds for (i, lmni) in enumerate(lmn_bp)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _dummy!(is,js,aijs,i,j)
            # _induction_sSS!(is,js,aijs,i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, sb; thresh, external)
            aij = _induction_sSS(T.(lmnj), T.(lmnb0), T.(lmni), r, wr, su, smfb0, sb; external)
            appendit!(is, js, aijs, i, j, aij; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            # _dummy!(is,js,aijs,i,j+np)
            # _induction_tSS!(is,js,aijs,i,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, s_mf; thresh, external)
            aij = _induction_tSS(T.(lmnj), T.(lmnb0), T.(lmni), r, wr, tu, smfb0, sb; external)
            appendit!(is, js, aijs, i, j + np, aij; thresh)
        end
    end

    @inbounds for (i, lmni) in enumerate(lmn_bt)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sST!(is, js, aijs, i + npb, j, T.(lmnj), T.(lmnb0), T.(lmni), r, wr, su, smfb0, tb; thresh)
            # _dummy!(is,js,aijs,i+npb,j)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj, mj, nj = lmnj
            conditions && nconditions && !ncondition(lb0, ni, nb0, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_tST!(is, js, aijs, i + npb, j + np, T.(lmnj), T.(lmnb0), T.(lmni), r, wr, tu, smfb0, tb; thresh)
            # _dummy!(is,js,aijs,i+npb,j+np)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatb, nmatu)
end

function rhs_induction_btor(N, m, lmnb0; ns=false, thresh::T=sqrt(eps()), conditions=true) where {T}
    su = s_in
    tu = t_in
    lb0, mb0, nb0 = lmnb0
    lmn_p = InviscidBasis.lmn_upol(N, m, ns)
    lmn_t = InviscidBasis.lmn_utor(N, m, ns)

    np = length(lmn_p)

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    npb = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]

    # r, wr = rquad(N+lb0+nb0+5)
    rwrs = [rquad(n + lb0 + nb0 + 1) for n in 1:N]

    for (i, lmni) in enumerate(lmn_bp)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sTS!(is, js, aijs, i, j, lmnj, lmnb0, lmni, r, wr, su, t_mf, s_mf; thresh)
        end
    end

    for (i, lmni) in enumerate(lmn_bt)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_p)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition1(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sTT!(is, js, aijs, i + npb, j, lmnj, lmnb0, lmni, r, wr, su, t_mf, t_mf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_t)
            lj, mj, nj = lmnj
            conditions && !ncondition(lb0, ni, nb0 + 1, nj) && continue
            conditions && !condition2(li, lb0, lj, mi, mb0, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_tTT!(is, js, aijs, i + npb, j + np, lmnj, lmnb0, lmni, r, wr, tu, t_mf, t_mf; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)
    nmatu = length(lmn_p) + length(lmn_t)

    return sparse(is, js, aijs, nmatb, nmatu)
end

function rhs_induction_upol(N, m, lmnu0; ns=false, thresh::T=sqrt(eps()), su=s_chen, smf=s_mf, tmf=t_mf, condition=true) where {T}
    lu0, mu0, nu0 = lmnu0

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    np = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    rwrs = [rquad(n + lu0 + nu0 + 1) for n in 1:N]

    for (i, lmni) in enumerate(lmn_bp)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sSS!(is, js, aijs, i, j, lmnu0, lmnj, lmni, r, wr, su, smf, smf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            !condition2u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sTS!(is, js, aijs, i, j + np, lmnu0, lmnj, lmni, r, wr, su, tmf, smf; thresh)
        end
    end

    for (i, lmni) in enumerate(lmn_bt)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sST!(is, js, aijs, i + np, j, lmnu0, lmnj, lmni, r, wr, su, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_sTT!(is, js, aijs, i + np, j + np, lmnu0, lmnj, lmni, r, wr, su, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)

    return sparse(is, js, aijs, nmatb, nmatb)
end

function rhs_induction_utor(N, m, lmnu0; ns=false, thresh::T=sqrt(eps()), tu=t_chen, smf=s_mf, tmf=t_mf, condition=true) where {T}
    lu0, mu0, nu0 = lmnu0

    lmn_bp = InsulatingMFBasis.lmn_bpol(N, m, ns)
    lmn_bt = InsulatingMFBasis.lmn_btor(N, m, ns)

    np = length(lmn_bp)

    is, js, aijs = Int[], Int[], Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    rwrs = [rquad(n + lu0 + nu0 + 1) for n in 1:N]

    for (i, lmni) in enumerate(lmn_bp)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_tSS!(is, js, aijs, i, j, lmnu0, lmnj, lmni, r, wr, tu, smf, smf; thresh)
        end
    end

    for (i, lmni) in enumerate(lmn_bt)
        li, mi, ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition1u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_tST!(is, js, aijs, i + np, j, lmnu0, lmnj, lmni, r, wr, tu, smf, tmf; thresh)
        end
        for (j, lmnj) in enumerate(lmn_bt)
            lj, mj, nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
            condition && !condition2u(lu0, li, lj, mu0, mi, mj) && continue
            r, wr = rwrs[min(N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            _induction_tTT!(is, js, aijs, i + np, j + np, lmnu0, lmnj, lmni, r, wr, tu, tmf, tmf; thresh)
        end
    end
    nmatb = length(lmn_bp) + length(lmn_bt)

    return sparse(is, js, aijs, nmatb, nmatb)
end


# function rhs_induction_btor_cond(N,m, lmnb0; ns = false,thresh::T = sqrt(eps())) where T
#     su = s_in
#     tu = t_in
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     tmf = t_in
#     smf = s_in
#     lmn_bp = InviscidBasis.lmn_upol(N,m,ns)
#     lmn_bt = InviscidBasis.lmn_utor(N,m,ns)

#     npb = length(lmn_bp)

#     is,js,aijs = Int[],Int[],Complex{T}[]

#     r, wr = rquad(N+lb0+nb0+5)

#     for (i,lmni) in enumerate(lmn_bp)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_sTS!(is,js,aijs,i,j,lmnj,lmnb0,lmni, r, wr, su, tmf, smf; thresh)
#         end
#     end

#     for (i,lmni) in enumerate(lmn_bt)
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_sTT!(is,js,aijs,i+npb,j,lmnj,lmnb0,lmni, r, wr, su, tmf, tmf; thresh)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             _induction_tTT!(is,js,aijs,i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, tmf, tmf; thresh)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(is,js,aijs,nmatb, nmatu)
# end


# function rhs_induction_bpol_dist(N,m, lmnb0; ns = false,thresh::T = sqrt(eps()), smfb0::Sf = s_mf) where {T,Sf}
#     su = s_in
#     tu = t_in
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     lmn_bp = InsulatingMFBasis.lmn_bpol(N,m,ns)
#     lmn_bt = InsulatingMFBasis.lmn_btor(N,m,ns)

#     npb = length(lmn_bp)


#     nt = nprocs()
#     is,js,aijs = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])

#     # r, wr = rquad(N+lb0+nb0+5)
#     rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]
#     @sync @distributed for i in shuffle(eachindex(lmn_bp))
#         lmni = lmn_bp[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0+1,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             # _dummy!(is,js,aijs,i,j)
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, s_mf; thresh)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             # _dummy!(is,js,aijs,i,j+np)
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_tSS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, s_mf; thresh)
#         end
#     end

#     @sync @distributed for i in shuffle(eachindex(lmn_bt))
#         lmni = lmn_bt[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, su, smfb0, t_mf; thresh)
#             # _dummy!(is,js,aijs,i+npb,j)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_tST!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j+np,T.(lmnj),T.(lmnb0),T.(lmni), r, wr, tu, smfb0, t_mf; thresh)
#             # _dummy!(is,js,aijs,i+npb,j+np)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatb, nmatu)
# end

# function rhs_induction_btor_dist(N,m, lmnb0; ns = false,thresh::T = sqrt(eps())) where T
#     su = s_in
#     tu = t_in
#     lb0,mb0,nb0 = lmnb0
#     lmn_p = InviscidBasis.lmn_upol(N,m,ns)
#     lmn_t = InviscidBasis.lmn_utor(N,m,ns)

#     np = length(lmn_p)

#     lmn_bp = InsulatingMFBasis.lmn_bpol(N,m,ns)
#     lmn_bt = InsulatingMFBasis.lmn_btor(N,m,ns)

#     npb = length(lmn_bp)

#     nt = nprocs()
#     is,js,aijs = distribute([Int[] for _ in 1:nt]),distribute([Int[] for _ in 1:nt]),distribute([Complex{T}[] for _ in 1:nt])


#     # r, wr = rquad(N+lb0+nb0+5)
#     rwrs = [rquad(n+lb0+nb0+1) for n in 1:N]

#     @sync @distributed for i in shuffle(eachindex(lmn_bp))
#         lmni = lmn_bp[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sTS!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, s_mf; thresh)
#         end
#     end

#     @sync @distributed for i in shuffle(eachindex(lmn_bt))
#         lmni = lmn_bt[i]
#         li,mi,ni = lmni
#         for (j, lmnj) in enumerate(lmn_p)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition1(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_sTT!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j,lmnj,lmnb0,lmni, r, wr, su, t_mf, t_mf; thresh)
#         end
#         for (j, lmnj) in enumerate(lmn_t)
#             lj,mj,nj = lmnj
#             !ncondition(lb0,ni,nb0,nj) && continue
#             !condition2(li,lb0,lj,mi,mb0,mj) && continue
#             r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
#             _induction_tTT!(first(localpart(is)),first(localpart(js)),first(localpart(aijs)),i+npb,j+np,lmnj,lmnb0,lmni, r, wr, tu, t_mf, t_mf; thresh)
#         end
#     end
#     nmatb = length(lmn_bp)+length(lmn_bt)
#     nmatu = length(lmn_p)+length(lmn_t)

#     return sparse(vcat(is...),vcat(js...),vcat(aijs...),nmatb, nmatu)
# end
