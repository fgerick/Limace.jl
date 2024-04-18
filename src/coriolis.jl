"""
$(TYPEDSIGNATURES)

Equation (113) in Ivers & Phillips (2008).
"""
function C(l, m)
    return (l^2 - 1) * √((l^2 - m^2) / (4l^2 - 1))
end

##fallbacks

function _coriolis_tt(::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_tt(T, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

function _coriolis_ss(::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_tt(T, lmna, lmnb, r, wr)
    return m * Ω / p(l) * aij
end

function _coriolis_st(::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
end

function _coriolis_ts(::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
end


@inline function _coriolis_poloidal_new(b::TB; Ω::T=2.0) where {TB<:Basis, T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    r, wr = rquad(b.N + 5)
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)


    for l in 1:lpmax(b)
        for m in intersect(b.m,-l:l)
            for n in nrange_p(b,l), n2 in nrange_p(b,l2)
                aij = _coriolis_ss(TB, (l, m, n), (l, m, n2), r, wr; Ω)
                appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
            end
            for l2 in  ((l == 1) ? (2,) : ((l == lpmax(b)) ? (l-1,) : (l - 1,l + 1)))
                for n in nrange_p(b,l), n2 in nrange_t(b,l2)
                    aij = _coriolis_st(TB, (l, m, n), (l2, m, n2), r, wr; Ω)
                    appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_t[(l2,m,n2)] + np, aij)
                end
            end
        end
    end

    return is, js, aijs
end

@inline function _coriolis_poloidal(b::TB, np, lmnp_l, lmnt_l; Ω::T=2.0) where {TB<:Basis, T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    r, wr = rquad(b.N + 5)

    for (l, ilmn) in enumerate(lmnp_l)
        for (i, _, m, n) in ilmn
            for (j, l2, m2, n2) in lmnp_l[l]
                if (m == m2)
                    aij = _coriolis_ss(TB, (l, m, n), (l2, m2, n2), r, wr; Ω)
                    appendit!(is, js, aijs, i, j, aij)
                end
            end
            lrange = (l == 1) ? [2] : ((l == lpmax(b)) ? [l - 1] : [l - 1, l + 1])
            for lmn2 in view(lmnt_l, lrange)
                for (j, l2, m2, n2) in lmn2
                    aij = _coriolis_st(TB, (l, m, n), (l2, m2, n2), r, wr; Ω)
                    appendit!(is, js, aijs, i, j + np, aij)
                end
            end
        end
    end

    return is, js, aijs
end

@inline function _coriolis_toroidal(b::TB, np, lmn_p_l, lmn_t_l; Ω::T=2.0) where {TB<:Basis, T}

    is, js, aijs = Int[], Int[], Complex{T}[]
	r, wr = rquad(b.N + 5)

    for (l, ilmn) in enumerate(lmn_t_l)
        for (i, _, m, n) in ilmn
            # l,m,n = T.((l,m,n))
            for (j, l2, m2, n2) in lmn_t_l[l]
                if (m == m2)
                    aij = _coriolis_tt(TB, (l,m,n), (l2,m2,n2), r, wr; Ω)
					appendit!(is, js, aijs, i + np, j + np, aij)
                end
            end
            lrange = (l == 1) ? [2] : ((l >= N - 1) ? [l - 1] : [l - 1, l + 1])
            for lmn2 in view(lmn_p_l, lrange)
                for (j, l2, m2, n2) in lmn2
                    aij = _coriolis_ts(TB, (l,m,n), (l2,m2,n2), r, wr; Ω)
					appendit!(is, js, aijs, i + np, j, aij)
                end
            end
        end
    end

    return is, js, aijs
end

function coriolis(b::TB; Ω::T=2.0) where {TB<:Basis,T}
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)
    lmnp_l = lmn_p_l(b)
    lmnt_l = lmn_t_l(b)

    np = length(lmnp)
    nt = length(lmnt)
    nu = np + nt

    is, js, aijs = _coriolis_poloidal(b, np, lmnp_l, lmnt_l; Ω)
    is2, js2, aijs2 = _coriolis_toroidal(b, np, lmnp_l, lmnt_l; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end

function coriolis_new(b::TB; Ω::T=2.0) where {TB<:Basis,T}
    lmnp = lmn_p(b)
    lmnt = lmn_t(b)
    lmnp_l = lmn_p_l(b)
    lmnt_l = lmn_t_l(b)

    np = length(lmnp)
    nt = length(lmnt)
    nu = np + nt

    is, js, aijs = _coriolis_poloidal_new(b; Ω)
    is2, js2, aijs2 = _coriolis_toroidal(b, np, lmnp_l, lmnt_l; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end