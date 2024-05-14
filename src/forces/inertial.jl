"""
$(TYPEDSIGNATURES)

"""
function _inertial_ss(b::T, lmna, lmnb, r,wr; external=false) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(T,b.V,la,ma,na,r)
    @inline _sb = r->s(T,b.V,lb,mb,nb,r)

    @inline f = r-> inners(_sa,_sb, la, r)

    aij = ∫dr(f,r,wr)

    if external && (la==lb) && (ma==mb)
        aij += _sa(b.V.r1)*_sb(b.V.r1)*la^2*(la+1)
    end
    return aij
end

"""
$(TYPEDSIGNATURES)

"""
function _inertial_tt(b::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(T,b.V,la,ma,na,r)
    @inline _tb = r->t(T,b.V,lb,mb,nb,r)

    @inline f = r-> innert(_ta,_tb, la, r)

    aij = ∫dr(f,r,wr)
    return aij
end

"""
$(TYPEDSIGNATURES)

"""
@inline function inertial(b::TB, ::Type{T}=Float64; external=false) where {TB<:Basis,T<:Number}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    nu = length(b)

    #m == m2 and only l==l2 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            for n in nrange_p_bc(b,l), n2 in nrange_p(b,l)
                aij = _inertial_ss(b, (l,m,n), (l,m,n2), r,wr; external)
                appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
            end
        end
    end

    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
                aij = _inertial_tt(b, (l,m,n), (l,m,n2), r,wr)
                appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
            end
        end
    end


    return sparse(is, js, aijs, nu, nu)
end

