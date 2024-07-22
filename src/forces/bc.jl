"""
$(TYPEDSIGNATURES)

"""
@inline function boundarycondition(b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    nbasis = length(b)

    #m == m2 and only l==l2 needs to be considered.
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

"""
$(TYPEDSIGNATURES)

"""
@inline function boundarycondition!(A, b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    nbasis = length(b)

    #m == m2 and only l==l2 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC()) 
                bcf = bcs_p(b) #tuple of evaluation functions
                npmax = last(nrange_p(b,l))
                for (i,f) in enumerate(bcf), n2 in nrange_p(b,l)
                    bij = f(l,n2)
                    A[lmn2k_p[(l,m,npmax-i+1)], lmn2k_p[(l,m,n2)]] = bij
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
                    A[lmn2k_t[(l,m,ntmax-i+1)] + _np, lmn2k_t[(l,m,n2)] + _np] = bij
                end
            end
        end
    end


    return nothing
end

"""
$(TYPEDSIGNATURES)

"""
@inline function zero_boundarycondition!(A, b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)

    #m == m2 and only l==l2 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC()) 
                bcf = bcs_p(b) #tuple of evaluation functions
                npmax = last(nrange_p(b,l))
                for (i,f) in enumerate(bcf), n2 in nrange_p(b,l)
                    A[lmn2k_p[(l,m,npmax-i+1)], lmn2k_p[(l,m,n2)]] = 0
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
                    A[lmn2k_t[(l,m,ntmax-i+1)] + _np, lmn2k_t[(l,m,n2)] + _np] = 0
                end
            end
        end
    end


    return nothing
end

@inline function zero_boundarycondition_row!(A, b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)

    #m == m2 and only l==l2 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC()) 
                bcf = bcs_p(b) #tuple of evaluation functions
                npmax = last(nrange_p(b,l))
                for (i,f) in enumerate(bcf)
                    A[lmn2k_p[(l,m,npmax-i+1)], :] .= 0
                end
            end
        end
    end

    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            if (b.BC == NoBC())
                bcf = bcs_t(b) #tuple of evaluation functions
                ntmax = last(nrange_t(b,l))
                for (i,f) in enumerate(bcf)
                    A[lmn2k_t[(l,m,ntmax-i+1)] + _np, :] .= 0
                end
            end
        end
    end


    return nothing
end