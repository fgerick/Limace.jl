const Tunion = Union{Insulating, Inviscid, Viscous, ThinWall} #all bases that can be used in this specialization

# only assemble cross terms that are neigboring up to max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))

@inline function _crossterm!(bbi::Basis{Ti,V}, U0::BasisElement{Basis{T0,V},PT,T}, bbj::Basis{Tj,V}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = U0.lmn
    for ni in nrangefi(bbi, li)
        njs_all = nrangefj(bbj,lj)
        for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(max(bbi.N,bbj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(Basis{T0,V}, Basis{Tj,V}, Basis{Ti,V}, bbi.V, U0.lmn, lmnj, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*U0.factor)
        end
    end
    return nothing
end

@inline function _crossterm!(bbi::Basis{Ti,V}, buj::Basis{Tj,V}, B0::BasisElement{Basis{T0,V},PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_uj, nrangefi, nrangefj, indf, EA; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bbi, li)
        # for nj in nrangefj(bbj, lj)
        njs_all = nrangefj(buj,lj)
        for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
            r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(Basis{Tj,V}, Basis{T0,V}, Basis{Ti,V}, bbi.V, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_uj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end


@inline function _crossterm_m_adamgaunt!(bbi::Basis{Ti,V}, buj::Basis{Tj,V}, B0::BasisElement{Basis{T0,V},PT,T}, is, js, aijs, i0, j0,
    li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = B0.lmn
    njs_all = nrangefj(buj,lj)
    As = ComplexF64[]
    for mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljsa
            push!(As,adamgaunt(lj,l0,li, mj, m0, mi))
        end
    end
    for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
        r, wr = rwrs[min(max(bbi.N,buj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
        aij = indf(Basis{Tj,V}, Basis{T0,V}, Basis{Ti,V}, bbi.V, (lj, 0, nj), B0.lmn,  (li, 0, ni), r, wr; kwargs...)
        iA = 1
        for mi in intersect(bbi.m, -li:li)
            mj = adamgaunt_mjs(mi, m0)
            ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
            if lj ∈ ljsa
                A = As[iA]
                iA +=1
                lmni = (li, mi, ni)
                lmnj = (lj, mj, nj)
                appendit!(is, js, aijs, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*A*B0.factor)
            end
        end
    end
    return nothing
end

@inline function _crossterm_m_elsasser!(bbi::Basis{Ti,V}, buj::Basis{Tj,V}, B0::BasisElement{Basis{T0,V},PT,T}, is, js, aijs, i0, j0,
    li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = B0.lmn
    njs_all = nrangefj(buj,lj)
    Es = ComplexF64[] 
    for mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljse
            push!(Es,elsasser(lj,l0,li, mj, m0, mi))
        end
    end
    for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
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
                appendit!(is, js, aijs, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*E*B0.factor)
            end
        end
    end
    return nothing
end



@inline function _crossterm_m_adamgaunt!(bbi::Basis{Ti,V}, U0::BasisElement{Basis{T0,V},PT,T}, buj::Basis{Tj,V}, is, js, aijs, i0, j0,
    li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = U0.lmn
    njs_all = nrangefj(buj,lj)
    As = ComplexF64[]
    for mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        ljsa = adamgaunt_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljsa
            push!(As,adamgaunt(l0,lj,li, m0, mj, mi))
        end
    end
    for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
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
                appendit!(is, js, aijs, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*A*U0.factor)
            end
        end
    end
    return nothing
end

@inline function _crossterm_m_elsasser!(bbi::Basis{Ti,V}, U0::BasisElement{Basis{T0,V},PT,T}, buj::Basis{Tj,V}, is, js, aijs, i0, j0,
    li, ni, lj, rwrs, lmn2k_bi, lmn2k_uj, nrangefj, lptmax, indf; kwargs...) where {Ti<:Tunion, Tj<:Tunion, T0<:Tunion, PT<:Helmholtz,T, V<:Volume}
    l0,m0,n0 = U0.lmn
    njs_all = nrangefj(buj,lj)
    Es = ComplexF64[] 
    for mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        ljse = elsasser_ljs(li, l0, mj, lptmax(buj))
        if lj ∈ ljse
            push!(Es,elsasser(l0,lj,li, m0, mj, mi))
        end
    end
    for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all)) #bandedness
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
                appendit!(is, js, aijs, lmn2k_bi[lmni]+i0, lmn2k_uj[lmnj]+j0, aij*E*U0.factor)
            end
        end
    end
    return nothing
end
