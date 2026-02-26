
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
@inline function _crossterm!(bi::TI, B0::BasisElement{T0,PT,T}, bj::TJ, is, js, aijs, lck, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(max(bi.N,bj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(T0, TJ, TI, bi.V, B0.lmn, lmnj, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lck, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fallback functions for `_crossterm!` term for `B0`. Write specialized function to include e.g. bandedness in `n`.
"""
@inline function _crossterm!(bi::TI, bj::TJ, B0::BasisElement{T0,PT,T}, is, js, aijs, lck, i0, j0,
    						 li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(max(bi.N,bj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(TJ, T0, TI, bi.V, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lck, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end

# """
# $(TYPEDSIGNATURES)

# """
# @inline function _crossterm!(bi::TI, B1::BasisElement{T1,PT,T}, B2::BasisElement{T2,PT,T}, is, aijs, lck, i0, 
# 							 li, mi, rwrs, lmn2k_bi, nrangefi, indf, EA; kwargs...) where {TI<:Basis,T1<:Basis,T2<:Basis,PT<:Helmholtz,T}
#     l1,m1,n1 = B1.lmn
#     l2,m2,n2 = B2.lmn
#     for ni in nrangefi(bi, li)
# 		r, wr = rwrs[min(bi.N, li ÷ 2 + ni + l1 + n1 + l2 + n2 + 1)]
# 		lmni = (li, mi, ni)
# 		aij = indf(T1, T2, TI, bi.V, B1.lmn, B2.lmn, lmni, r, wr; kwargs...)*EA
# 		appendit!(is, aijs, lck, lmn2k_bi[lmni] + i0, aij*B1.factor*B2.factor)
#     end
#     return nothing
# end

@inline function _crossterm_m_adamgaunt!(bbi::Basis{Ti,V}, buj::Basis{Tj,V}, B0::BasisElement{Basis{T0,V},PT,T}, is, js, aijs, lck, i0, j0,
    									 li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti, Tj, T0, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = B0.lmn
    As = ComplexF64[]
    for mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljsa
            push!(As,adamgaunt(lj,l0,li, mj, m0, mi))
        end
    end
    for nj in nrangefj(buj,lj)
        r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
        aij = indf(Basis{Tj,V}, Basis{T0,V}, Basis{Ti,V}, bbi.V, (lj, 0, nj), B0.lmn,  (li, 0, ni), r, wr; kwargs...)
        iA = 1
        for mi in intersect(bbi.m, -li:li)
            mj = adamgaunt_mjs(mi, m0)
            ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
            if lj ∈ ljsa
                A = As[iA]
                iA+=1
                lmni = (li, mi, ni)
                lmnj = (lj, mj, nj)
                appendit!(is, js, aijs, lck, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*A*B0.factor)
            end
        end
    end
    return nothing
end

@inline function _crossterm_m_elsasser!(bbi::Basis{Ti,V}, buj::Basis{Tj,V}, B0::BasisElement{Basis{T0,V},PT,T}, is, js, aijs, lck, i0, j0,
    									li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj,lptmax, indf; kwargs...) where {Ti, Tj, T0, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = B0.lmn
    Es = ComplexF64[]
    for mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljse
            push!(Es,elsasser(lj,l0,li, mj, m0, mi))
        end
    end
    for nj in nrangefj(buj,lj)
        r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
        aij = indf(Basis{Tj,V}, Basis{T0,V}, Basis{Ti,V}, bbi.V, (lj, 0, nj), B0.lmn,  (li, 0, ni), r, wr; kwargs...)
        iE = 1
        for mi in intersect(bbi.m, -li:li)
            mj = elsasser_mjs(mi, m0)
            ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
            if lj ∈ ljse
                E = Es[iE]
                iE+=1
                lmni = (li, mi, ni)
                lmnj = (lj, mj, nj)
                appendit!(is, js, aijs, lck, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*E*B0.factor)
            end
        end
    end
    return nothing
end


@inline function _crossterm_m_adamgaunt!(bbi::Basis{Ti,V}, U0::BasisElement{Basis{T0,V},PT,T}, buj::Basis{Tj,V}, is, js, aijs, lck, i0, j0,
    									 li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti, Tj, T0, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = U0.lmn
    As = ComplexF64[]
    for mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljsa
            push!(As,adamgaunt(l0,lj,li, m0, mj, mi))
        end
    end
    for nj in nrangefj(buj,lj)
        r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
        aij = indf(Basis{T0,V}, Basis{Tj,V}, Basis{Ti,V}, bbi.V, U0.lmn, (lj, 0, nj), (li, 0, ni), r, wr; kwargs...)
        iA = 1
        for mi in intersect(bbi.m, -li:li)
            mj = adamgaunt_mjs(mi, m0)
            ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
            if lj ∈ ljsa
                A = As[iA]
                iA +=1
                lmni = (li, mi, ni)
                lmnj = (lj, mj, nj)
                appendit!(is, js, aijs, lck, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*A*U0.factor)
            end
        end
    end
    return nothing
end

@inline function _crossterm_m_elsasser!(bbi::Basis{Ti,V}, U0::BasisElement{Basis{T0,V},PT,T}, buj::Basis{Tj,V}, is, js, aijs, lck, i0, j0,
    li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti, Tj, T0, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = U0.lmn
    Es = ComplexF64[] 
    for mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljse
            push!(Es,elsasser(l0,lj,li, m0, mj, mi))
        end
    end
    for nj in nrangefj(buj,lj)
        r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
        aij = indf(Basis{T0,V}, Basis{Tj,V}, Basis{Ti,V}, bbi.V, U0.lmn,  (lj, 0, nj), (li, 0, ni), r, wr; kwargs...)
        iE = 1
        for mi in intersect(bbi.m, -li:li)
            mj = elsasser_mjs(mi, m0)
            ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
            if lj ∈ ljse
                E = Es[iE]
                iE+=1
                lmni = (li, mi, ni)
                lmnj = (lj, mj, nj)
                appendit!(is, js, aijs, lck, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*E*U0.factor)
            end
        end
    end
    return nothing
end

