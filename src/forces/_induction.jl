
#∫Sᵢ⋅∇×(pⱼ×S₀) dV
function _induction_poloidal_poloidal_poloidal(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_p_uj = lmn2k_p_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]
    # r, wr = rquad(N + l0+n0+1)

    nbi = length(lmn2k_p_bi)
    nuj = length(lmn2k_p_uj)

    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj)
            _induction_poloidal_poloidal_poloidal!(bbi, buj, B0, is, js, aijs, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj; external)
        end
    end

    return sparse(is, js, aijs, nbi, nuj)
end

# function _induction_poloidal_poloidal_poloidal!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs,
#     li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
#     for ni in nrange_p(bbi, li)
#         for nj in nrange_p(buj, lj)
#             # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
#             r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
#             lmni = (li, mi, ni)
#             lmnj = (lj, mj, nj)
#             aij = _induction_sSS(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr; external)
#             appendit!(is, js, aijs, lmn2k_p_bi[lmni], lmn2k_p_uj[lmnj], aij)
#         end
#     end
#     return nothing
# end

#∫Sᵢ⋅∇×(p₀×Sⱼ) dV
function _induction_poloidal_poloidal_poloidal(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_p_bj = lmn2k_p_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]
    # r, wr = rquad(N + l0+n0+1)

    nbi = length(lmn2k_p_bi)
    nbj = length(lmn2k_p_bj)

    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, -m0)
        for lj in adamgaunt_ljs(li, l0, mj)
            _induction_poloidal_poloidal_poloidal!(bbi, U0, bbj, is, js, aijs, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj; external)
        end
    end

    return sparse(is, js, aijs, nbi, nbj)
end

# function _induction_poloidal_poloidal_poloidal!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs,
#     li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
#     for ni in nrange_p(bbi, li)
#         for nj in nrange_p(bbj, lj)
#             # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
#             r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
#             lmni = (li, mi, ni)
#             lmnj = (lj, mj, nj)
#             aij = _induction_sSS(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr; external)
#             appendit!(is, js, aijs, lmn2k_p_bi[lmni], lmn2k_p_bj[lmnj], aij)
#         end
#     end
#     return nothing
# end


#∫Sᵢ⋅∇×(pⱼ×T₀) dV
function _induction_poloidal_poloidal_toroidal(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_p_uj = lmn2k_p_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]
    # r, wr = rquad(N + l0+n0+1)

    nbi = length(lmn2k_p_bi)
    nuj = length(lmn2k_p_uj)

    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj)
            _induction_poloidal_poloidal_toroidal!(bbi, buj, B0, is, js, aijs, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj)
        end
    end

    return sparse(is, js, aijs, nbi, nuj)
end

# function _induction_poloidal_poloidal_toroidal!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs,
#     li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
#     for ni in nrange_p(bbi, li)
#         for nj in nrange_p(buj, lj)
#             # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
#             r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
#             lmni = (li, mi, ni)
#             lmnj = (lj, mj, nj)
#             aij = _induction_sTS(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr)
#             appendit!(is, js, aijs, lmn2k_p_bi[lmni], lmn2k_p_uj[lmnj], aij)
#         end
#     end
#     return nothing
# end


#∫Sᵢ⋅∇×(p₀×Tⱼ) dV
function _induction_poloidal_poloidal_toroidal(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    is, js, aijs = Int[], Int[], Complex{T}[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1) for n in 1:N]
    # r, wr = rquad(N + l0+n0+1)

    nbi = length(lmn2k_p_bi)
    nbj = length(lmn2k_p_bj)

    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, -m0)
        for lj in elasser_ljs(li, l0, mj)
            _induction_poloidal_poloidal_toroidal!(bbi, U0, bbj, is, js, aijs, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj; external)
        end
    end

    return sparse(is, js, aijs, nbi, nbj)
end

# function _induction_poloidal_poloidal_toroidal!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs,
#     li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
#     for ni in nrange_p(bbi, li)
#         for nj in nrange_p(bbj, lj)
#             # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
#             r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
#             lmni = (li, mi, ni)
#             lmnj = (lj, mj, nj)
#             aij = _induction_sTS(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr)
#             appendit!(is, js, aijs, lmn2k_p_bi[lmni], lmn2k_p_bj[lmnj], aij)
#         end
#     end
#     return nothing
# end











#∫Tᵢ⋅∇×(pⱼ×S₀) dV
# function _induction_toroidal_poloidal_poloidal(bbi::Basis, buj::Basis, B0::BasisElement)
# end

function _induction_toroidal_poloidal_poloidal!(bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_p(buj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_sST(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_p_uj[lmnj], aij)
        end
    end
    return nothing
end


#∫Tᵢ⋅∇×(p₀×Sⱼ) dV
# function _induction_toroidal_poloidal_poloidal(bbi::Basis, U0::BasisElement, bbj::Basis)
# end

function _induction_toroidal_poloidal_poloidal!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_p(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_sST(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_p_bj[lmnj], aij)
        end
    end
    return nothing
end

#∫Tᵢ⋅∇×(pⱼ×T₀) dV
# function _induction_toroidal_poloidal_toroidal(bbi::Basis, buj::Basis, B0::BasisElement)
# end

function _induction_toroidal_poloidal_toroidal!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_p(buj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_sTT(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_p_uj[lmnj], aij)
        end
    end
    return nothing
end

#∫Tᵢ⋅∇×(p₀×Tⱼ) dV
# function _induction_toroidal_poloidal_toroidal(bbi::Basis, U0::BasisElement, bbj::Basis)
# end

function _induction_toroidal_poloidal_toroidal!(bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_t(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_sTT(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_t_bj[lmnj], aij)
        end
    end
    return nothing
end


#∫Tᵢ⋅∇×(qⱼ×T₀) dV
# function _induction_toroidal_toroidal_toroidal(bbi::Basis, buj::Basis, B0::BasisElement)
# end

function _induction_toroidal_toroidal_toroidal!(bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_t(buj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_tTT(TJ, T0, TI, lmnj, B0.lmn, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_t_uj[lmnj], aij)
        end
    end
    return nothing
end

#∫Tᵢ⋅∇×(q₀×Tⱼ) dV
# function _induction_toroidal_toroidal_toroidal(bbi::Basis, U0::BasisElement, bbj::Basis)
# end

function _induction_toroidal_toroidal_toroidal!(bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ, is, js, aijs,
    li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    for ni in nrange_t(bbi, li)
        for nj in nrange_t(bbj, lj)
            # for nj in max(ni-(n0+1)-l0,first(njs_all)):min(ni+n0+1+l0,last(njs_all))
            r, wr = rwrs[min(bbi.N, li ÷ 2 + ni + lj ÷ 2 + nj + 1)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = _induction_tTT(T0, TJ, TI, U0.lmn, lmnj, lmni, r, wr)
            appendit!(is, js, aijs, lmn2k_t_bi[lmni], lmn2k_t_bj[lmnj], aij)
        end
    end
    return nothing
end


