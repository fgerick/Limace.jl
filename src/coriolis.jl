"""
$(TYPEDSIGNATURES)

Equation (113) in Ivers & Phillips (2008).
"""
function C(l, m)
    return (l^2 - 1) * √((l^2 - m^2) / (4l^2 - 1))
end

# Ivers & Phillips (2008) eq. (25)
function _∂ll(f,l,l1,r) 
    @assert l1 ∈ (l-1, l+1)
    if l1 == l-1
        return ∂(f,r) + (l+1)/r*f(r) 
    elseif l1 == l+1
        return ∂(f,r)-l/r*f(r)
    end
end

##fallbacks

function _coriolis_tt(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_tt(b, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

function _coriolis_ss(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_ss(b, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

function _coriolis_st(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    @inline _sa = r->s(T,la,ma,na,r)
    @inline _tb = r->t(T,lb,mb,nb,r)

    if lb == la-1

        @inline f1 = r-> C(la,ma)*_∂ll(_tb,lb,la,r)
        @inline f = r -> innert(_sa,f1, la, r)
        aij = ∫dr(f,r,wr)

        return Ω / p(la) * aij
    elseif lb == la+1

        @inline f1 = r-> C(la+1,ma)*_∂ll(_tb,lb,la,r)
        @inline f = r-> innert(_sa,f1, la, r)
        aij = ∫dr(f,r,wr)

        return Ω / p(la) * aij
    else
        return nothing
    end
end

function _coriolis_ts(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    @inline _ta = r->t(T,la,ma,na,r)
    @inline _sb = r->s(T,lb,mb,nb,r)

    if lb == la-1

        @inline f1 = r-> C(la,ma)*_∂ll(_sb,lb,la,r)
        @inline f = r -> innert(_ta,f1, la, r)
        aij = ∫dr(f,r,wr)

        return Ω / p(la) * aij
    elseif lb == la+1

        @inline f1 = r-> C(la+1,ma)*_∂ll(_sb,lb,la,r)
        @inline f = r-> innert(_ta,f1, la, r)
        aij = ∫dr(f,r,wr)

        return Ω / p(la) * aij
    else
        return nothing
    end
end

function _coriolis_poloidal_poloidal!(b::Basis, is, js, aijs, lmn2k_p, l, m, r, wr, Ω)
    for n in nrange_p(b, l), n2 in nrange_p(b, l)
        aij = _coriolis_ss(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_p[(l, m, n)], lmn2k_p[(l, m, n2)], aij)
    end
    return nothing
end

function _coriolis_poloidal_toroidal!(b::Basis, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)
    for n in nrange_p(b, l), n2 in nrange_t(b, l2)
        aij = _coriolis_st(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_p[(l, m, n)], lmn2k_t[(l2, m, n2)] + _np, aij)
    end
    return nothing
end


function _coriolis_toroidal_toroidal!(b::Basis, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω)
    for n in nrange_t(b, l), n2 in nrange_t(b, l)
        aij = _coriolis_tt(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_t[(l, m, n)] + _np, lmn2k_t[(l, m, n2)] + _np, aij)
    end
    return nothing
end

function _coriolis_toroidal_poloidal!(b::Basis, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)
    for n in nrange_t(b, l), n2 in nrange_p(b, l2)
        aij = _coriolis_ts(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_t[(l, m, n)] + _np, lmn2k_p[(l2, m, n2)], aij)
    end
    return nothing
end

@inline function _coriolis_poloidal(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5)

    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_poloidal_poloidal!(b, is, js, aijs, lmn2k_p, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l == ltmax(b)) ? (l - 1,) : (l - 1, l + 1)))
                if l2 >= abs(m)
                    _coriolis_poloidal_toroidal!(b, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end


@inline function _coriolis_toroidal(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5)

    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_toroidal_toroidal!(b, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l == lpmax(b)) ? (l - 1,) : (l - 1, l + 1)))
                if abs(m) <= l2
                    _coriolis_toroidal_poloidal!(b, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end

function coriolis(b::TB; Ω::T=2.0) where {TB<:Basis,T}
    nu = length(b)

    is, js, aijs = _coriolis_poloidal(b; Ω)
    is2, js2, aijs2 = _coriolis_toroidal(b; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end
