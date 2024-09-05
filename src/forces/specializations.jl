@inline function _crossterm!(bbi::Basis{Ti}, U0::BasisElement{Basis{T0},PT,T}, bbj::Basis{Tj}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {Ti<:Union{Insulating,Inviscid,Viscous,ThinWall}, Tj<:Union{Insulating,Inviscid,Viscous,ThinWall}, T0<:Union{Insulating,Inviscid,Viscous,ThinWall}, PT<:Helmholtz,T}
    l0,m0,n0 = U0.lmn
    for ni in nrangefi(bbi, li)
        njs_all = nrangefj(bbj,lj)
        for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(max(bbi.N,bbj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(Basis{T0}, Basis{Tj}, Basis{Ti}, bbi.V, U0.lmn, lmnj, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*U0.factor)
        end
    end
    return nothing
end

@inline function _crossterm!(bbi::Basis{Ti}, buj::Basis{Tj}, B0::BasisElement{Basis{T0},PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_uj, nrangefi, nrangefj, indf, EA; kwargs...) where {Ti<:Union{Insulating,Inviscid,Viscous,ThinWall}, Tj<:Union{Insulating,Inviscid,Viscous,ThinWall}, T0<:Union{Insulating,Inviscid,Viscous,ThinWall}, PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bbi, li)
        # for nj in nrangefj(bbj, lj)
        njs_all = nrangefj(buj,lj)
        for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
            r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(Basis{Tj}, Basis{T0}, Basis{Ti}, bbi.V, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_uj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end
