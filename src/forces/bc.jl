"""
$(TYPEDSIGNATURES)

"""
@inline function boundarycondition(b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    nbasis = length(b)

    #m == m2 and only l==l2 needs to be considered.
    # k=1
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC()) 
                bcf = bcs_p(b) #tuple of evaluation functions
                npmax = last(nrange_p(b,l))
                for (i,f) in enumerate(bcf), n2 in nrange_p(b,l)
                    bij = f(l,n2)
                    appendit!(is, js, aijs, lmn2k_p[(l,m,npmax-i+1)], lmn2k_p[(l,m,n2)], bij)
                end
            end
        end
    end

    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC())
                bcf = bcs_t(b) #tuple of evaluation functions
                ntmax = last(nrange_t(b,l))
                for (i,f) in enumerate(bcf), n2 in nrange_t(b,l)
                    bij = f(l,n2)
                    appendit!(is, js, aijs, lmn2k_t[(l,m,ntmax-i+1)] + _np, lmn2k_t[(l,m,n2)] + _np, bij)
                end
            end
        end
    end


    return sparse(is, js, aijs, nbasis, nbasis)
end