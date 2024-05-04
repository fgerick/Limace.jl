"""
$(TYPEDSIGNATURES)

"""
function _diffusion_ss(b::T, lmna, lmnb, r,wr; external=false) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(T,b.V,la,ma,na,r)
    @inline _sb = r->s(T,b.V,lb,mb,nb,r)

    if external
        @inline f1 = r -> D(_sa,la,r)
        @inline f2 = r -> D(_sb,lb,r)
        @inline f = r-> -innert(f1,f2, lb, r)
    else
        @inline f1 = r -> D(_sa,la,r)
        @inline f = r-> inners(_sb,f1, lb, r)
    end

    aij = ∫dr(f,r,wr)
    return aij
end

"""
$(TYPEDSIGNATURES)

"""
function _diffusion_tt(b::T, lmna, lmnb, r,wr; external = false) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(T,b.V,la,ma,na,r)
    @inline _tb = r->t(T,b.V,lb,mb,nb,r)

    if external
        @inline f = r-> -inners(_tb,_ta, lb, r)
    else
        @inline f1 = r -> D(_ta,la,r)
        @inline f = r-> innert(_tb,f1, lb, r)

    end

    aij = ∫dr(f,r,wr)

    return aij
end


"""
$(TYPEDSIGNATURES)

"""
@inline function diffusion(b::TB, ::Type{T}=Float64; applyBC = true, external=false) where {TB<:Basis,T<:Number}

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
            for n in nrange_p_bc(b,l), n2 in nrange_p(b,l)
                aij = _diffusion_ss(b, (l,m,n), (l,m,n2), r,wr; external)
                # k==1 && @show aij
                # k+=1
                appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
            end
            if (b.BC == NoBC()) && applyBC
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
            for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
                aij = _diffusion_tt(b, (l,m,n), (l,m,n2), r,wr; external)
                appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
            end
            if (b.BC == NoBC()) && applyBC
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
@inline function diffusion_threaded(b::TB, ::Type{T}=Float64) where {TB<:Basis,T<:Number}

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5,b.V)
    nbasis = length(b)

    #m == m2 and only l==l2 needs to be considered.
    @sync begin
        for l in 1:lpmax(b)
            for m in intersect(b.m, -l:l)
                Threads.@spawn begin
                    id = Threads.threadid()
                    for n in nrange_p_bc(b,l), n2 in nrange_p(b,l)
                        aij = _diffusion_ss(b, (l,m,n), (l,m,n2), r,wr)
                        appendit!(is[id], js[id], aijs[id], lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
                    end
                    if b.BC == NoBC()
                        bcf = bcs_p(b) #tuple of evaluation functions
                        npmax = last(nrange_p(b,l))
                        for (i,f) in enumerate(bcf), n2 in nrange_p(b,l)
                            bij = f(l,n2)
                            appendit!(is[id], js[id], aijs[id], lmn2k_p[(l,m,npmax)]-i+1, lmn2k_p[(l,m,n2)], bij)
                        end
                    end
                end
            end
        end

        for l in 1:ltmax(b)
            for m in intersect(b.m, -l:l)
                Threads.@spawn begin
                    id = Threads.threadid()
                    for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
                        aij = _diffusion_tt(b, (l,m,n), (l,m,n2), r,wr)
                        appendit!(is[id], js[id], aijs[id], lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
                    end
                    if b.BC == NoBC()
                        bcf = bcs_t(b) #tuple of evaluation functions
                        ntmax = last(nrange_t(b,l))
                        for (i,f) in enumerate(bcf), n2 in nrange_t(b,l)
                            bij = f(l,n2)
                            appendit!(is[id], js[id], aijs[id], lmn2k_t[(l,m,ntmax)]-i+1 + _np, lmn2k_t[(l,m,n2)] + _np, bij)
                        end
                    end
                end
            end
        end
    end


    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nbasis, nbasis)
end


