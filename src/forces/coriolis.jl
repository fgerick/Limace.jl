"""
$(TYPEDSIGNATURES)

Equation (113) in Ivers & Phillips (2008).
"""
function C(l, m)
    return (l^2 - 1) * √((l^2 - m^2) / (4l^2 - 1))
end

##fallbacks
"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between toroidal and toroidal component, explicitly calculating the quadrature.
"""
function _coriolis_tt(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_tt(b, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between poloidal and poloidal component, explicitly calculating the quadrature.
"""
function _coriolis_ss(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    l, m, n = lmna
    aij = _inertial_ss(b, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between poloidal and toroidal component, explicitly calculating the quadrature. 
Following eq. (114) in Ivers & Phillips (2008).
"""
function _coriolis_st(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    @inline _sa = r->s(T,b.V,la,ma,na,r)
    @inline _tb = r->t(T,b.V,lb,mb,nb,r)

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

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between toroidal and poloidal component, explicitly calculating the quadrature.
Following eq. (112) in Ivers & Phillips (2008).
"""
function _coriolis_ts(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    @inline _ta = r->t(T,b.V,la,ma,na,r)
    @inline _sb = r->s(T,b.V,lb,mb,nb,r)

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

function _coriolis_poloidal_poloidal!(b::T, is, js, aijs, lmn2k_p, l, m, r, wr, Ω; applyBC=true) where T<:Basis
    for n in nrange_p_bc(b, l), n2 in nrange_p(b, l)
        aij = _coriolis_ss(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_p[(l, m, n)], lmn2k_p[(l, m, n2)], aij)
    end

    if (b.BC == NoBC()) && applyBC
        bcf = bcs_p(b) #tuple of evaluation functions
        npmax = last(nrange_p(b,l))
        for (i,f) in enumerate(bcf), n2 in nrange_p(b,l)
            bij = f(l,n2)
            appendit!(is, js, aijs, lmn2k_p[(l,m,npmax-i+1)], lmn2k_p[(l,m,n2)], bij)
        end
    end
    return nothing
end

function _coriolis_poloidal_toroidal!(b::T, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω) where T<:Basis
    for n in nrange_p_bc(b, l), n2 in nrange_t(b, l2)
        aij = _coriolis_st(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_p[(l, m, n)], lmn2k_t[(l2, m, n2)] + _np, aij)
    end
    return nothing
end


function _coriolis_toroidal_toroidal!(b::T, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω; applyBC=true) where T<:Basis
    for n in nrange_t_bc(b, l), n2 in nrange_t(b, l)
        aij = _coriolis_tt(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_t[(l, m, n)] + _np, lmn2k_t[(l, m, n2)] + _np, aij)
    end
    if (b.BC == NoBC()) && applyBC
        bcf = bcs_t(b) #tuple of evaluation functions
        ntmax = last(nrange_t(b,l))
        for (i,f) in enumerate(bcf), n2 in nrange_t(b,l)
            bij = f(l,n2)
            appendit!(is, js, aijs, lmn2k_t[(l,m,ntmax-i+1)] + _np, lmn2k_t[(l,m,n2)] + _np, bij)
        end
    end
    return nothing
end

function _coriolis_toroidal_poloidal!(b::T, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω) where T<:Basis
    for n in nrange_t_bc(b, l), n2 in nrange_p(b, l2)
        aij = _coriolis_ts(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lmn2k_t[(l, m, n)] + _np, lmn2k_p[(l2, m, n2)], aij)
    end
    return nothing
end

@inline function _coriolis_poloidal(b::Basis; Ω::T=2.0, applyBC=true) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)


    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_poloidal_poloidal!(b, is, js, aijs, lmn2k_p, l, m, r, wr, Ω; applyBC)
            for l2 in ((l == 1) ? (2,) : ((l+1 > ltmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_poloidal_toroidal!(b, is, js, aijs, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end


@inline function _coriolis_toroidal(b::Basis; Ω::T=2.0, applyBC=true) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_toroidal_toroidal!(b, is, js, aijs, _np, lmn2k_t, l, m, r, wr, Ω; applyBC)
            for l2 in ((l == 1) ? (2,) : ((l+1 > lpmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_toroidal_poloidal!(b, is, js, aijs, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end

function coriolis(b::TB; Ω::T=2.0, applyBC=true) where {TB<:Basis,T}
    nu = length(b)

    is, js, aijs = _coriolis_poloidal(b; Ω, applyBC)
    is2, js2, aijs2 = _coriolis_toroidal(b; Ω, applyBC)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end
