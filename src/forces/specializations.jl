@inline function __induction!(bbi::Basis{Insulating}, bbj::Basis{Inviscid}, B0::BasisElement{T0,PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_uj, nrangefi, nrangefj, indf, EA; kwargs...) where {T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bbi, li)
        # for nj in nrangefj(bbj, lj)
        njs_all = nrangefj(bbj,lj)
        for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(Basis{Inviscid}, T0, Basis{Insulating}, lmnj, B0.lmn, lmni, r, wr; kwargs...)
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_uj[lmnj] + j0, aij*EA)
        end
    end
    return nothing
end
