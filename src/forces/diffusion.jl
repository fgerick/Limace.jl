"""
$(TYPEDSIGNATURES)

Diffusion term of the poloidal components. Set `external=true` to include the contribution 
of an external poloidal magnetic field matching at the surface.
"""
function _diffusion_ss(b::T, lmna, lmnb, r,wr; external=false) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(T,b.V,la,ma,na,r)
    @inline _sb = r->s(T,b.V,lb,mb,nb,r)

    @inline f1 = r -> D(_sb,lb,r)
    @inline f = r-> inners(_sa,f1, la, r)

    aij = ∫dr(f,r,wr)


    if external
        r1 = b.V.r1
        aij += la^2*(la+1)*(∂(r->∂(r->r*_sb(r),r),r1)-lb*(lb+1)/r1*_sb(r1))*_sa(r1)
    end
    
    return aij
end

"""
$(TYPEDSIGNATURES)

Diffusion term of the toroidal components.
"""
function _diffusion_tt(b::T, lmna, lmnb, r,wr) where T<:Basis
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(T,b.V,la,ma,na,r)
    @inline _tb = r->t(T,b.V,lb,mb,nb,r)

    @inline f1 = r -> D(_tb,lb,r)
    @inline f = r-> innert(_ta,f1, la, r)

    aij = ∫dr(f,r,wr)

    return aij
end


"""
$(TYPEDSIGNATURES)

Compute the Galerkin projection matrix of the basis `b` onto the vector Laplacian. When keyword `external=true`, 
the integral is computed over all space, assuming continuity of the poloidal field and a scalar potential in the exterior domain.
"""
@inline function diffusion(b::Basis; external=false)
    T = typeof(b.V.r1)
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
                appendit!(is, js, aijs, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
            end
        end
    end

    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
                aij = _diffusion_tt(b, (l,m,n), (l,m,n2), r,wr)
                appendit!(is, js, aijs, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
            end
        end
    end


    return sparse(is, js, aijs, nbasis, nbasis)
end

"""
$(TYPEDSIGNATURES)

"""
@inline function diffusion_threaded(b::Basis; external=false)
    T = typeof(b.V.r1)

    is, js, aijs = Int[], Int[], Complex{T}[]
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [Complex{T}[] for _ in 1:_nt]

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
                        aij = _diffusion_ss(b, (l,m,n), (l,m,n2), r,wr; external)
                        appendit!(is[id], js[id], aijs[id], lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
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
                end
            end
        end
    end


    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nbasis, nbasis)
end


