"""
$(TYPEDEF)

- `basis::Basis{TB}`: The basis for the diffusion operator.
- `factor::T`: : A scalar factor that multiplies the diffusion operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the diffusion operator.
- `preassembled::Bool`: A flag indicating whether the diffusion operator has been preassembled.

"""
mutable struct Diffusion{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Diffusion(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Diffusion(b, ComplexF64(factor), mat, false)
end

function assemble!(f::Diffusion; kwargs...)
    f.mat = diffusion(f.basis; kwargs...)
    f.preassembled = true
    return f.mat
end

"""
$(TYPEDSIGNATURES)

Diffusion term of the poloidal components. Set `external=true` to include the contribution 
of an external poloidal magnetic field matching at the surface.
"""
function _diffusion_ss(bi::Ti, bj::Tj, lmna, lmnb, r,wr; external=false) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(Ti,bi.V,la,ma,na,r)
    @inline _sb = r->s(Tj,bj.V,lb,mb,nb,r)

    @inline f1 = r -> D(_sb,lb,r)
    @inline f = r-> inners(_sa,f1, la, r)

    aij = ∫dr(f,r,wr)


    if external
        r1 = bi.V.r1
        __sb, _dsb, _d2sb = derivatives012(_sb,r1)
        aij += la^2*(la+1)*((r1*_d2sb+2*_dsb)-lb*(lb+1)/r1*__sb)*_sa(r1)
    end
    
    return aij
end

"""
$(TYPEDSIGNATURES)

Diffusion term of the toroidal components.
"""
function _diffusion_tt(bi::Ti, bj::Tj, lmna, lmnb, r,wr) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(Ti,bi.V,la,ma,na,r)
    @inline _tb = r->t(Tj,bj.V,lb,mb,nb,r)

    @inline f1 = r -> D(_tb,lb,r)
    @inline f = r-> innert(_ta,f1, la, r)

    aij = ∫dr(f,r,wr)

    return aij
end


@inline function _diffusion(::Val{false}, bi::Basis, bj::Basis; external=false)
    T = promote_type(typeof(bi.V.r1),typeof(bj.V.r1))
    is, js, aijs = Int[], Int[], Complex{T}[]
    lck = ReentrantLock()
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    nbasisi = length(bi)
    nbasisj = length(bj)

    #m == m2 and only l==l2 needs to be considered.
    # k=1
    for l in 1:lpmax(bi)
        for m in intersect(bi.m, -l:l)
            for n in nrange_p_bc(bi,l), n2 in nrange_p(bj,l)
                aij = _diffusion_ss(bi,bj, (l,m,n), (l,m,n2), r,wr; external)
                appendit!(is, js, aijs, lck, lmn2k_pi[(l,m,n)], lmn2k_pj[(l,m,n2)], aij)
            end
        end
    end

    for l in 1:ltmax(bi)
        for m in intersect(bi.m, -l:l)
            for n in nrange_t_bc(bi,l), n2 in nrange_t(bj,l)
                aij = _diffusion_tt(bi,bj, (l,m,n), (l,m,n2), r,wr)
                appendit!(is, js, aijs, lck, lmn2k_ti[(l,m,n)] + _npi, lmn2k_tj[(l,m,n2)] + _npj, aij)
            end
        end
    end


    return sparse(is, js, aijs, nbasisi, nbasisj)
end

_diffusion(::Val{false}, b::Basis; external=false) = _diffusion(Val(false), b, b; external)

@inline function _diffusion(::Val{true}, bi::Basis, bj::Basis; external=false)
    T = promote_type(typeof(bi.V.r1),typeof(bj.V.r1))
    is, js, aijs = Int[], Int[], Complex{T}[]
    lck = ReentrantLock()
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    nbasisi = length(bi)
    nbasisj = length(bj)

    #m == m2 and only l==l2 needs to be considered.
    @sync begin
        for l in 1:lpmax(bi)
            for m in intersect(bi.m, -l:l)
                Threads.@spawn begin
                    for n in nrange_p_bc(bi,l), n2 in nrange_p(bj,l)
                        aij = _diffusion_ss(bi, bj, (l,m,n), (l,m,n2), r,wr; external)
                        appendit!(is, js, aijs, lck, lmn2k_pi[(l,m,n)], lmn2k_pj[(l,m,n2)], aij)
                    end
                end
            end
        end

        for l in 1:ltmax(bi)
            for m in intersect(bi.m, -l:l)
                Threads.@spawn begin
                    for n in nrange_t_bc(bi,l), n2 in nrange_t(bj,l)
                        aij = _diffusion_tt(bi, bj, (l,m,n), (l,m,n2), r,wr)
                        appendit!(is, js, aijs, lck, lmn2k_ti[(l,m,n)] + _npi, lmn2k_tj[(l,m,n2)] + _npj, aij)
                    end
                end
            end
        end
    end


    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nbasisi, nbasisj)
end

_diffusion(::Val{true}, b::Basis; external=false) = _diffusion(Val(true), b, b; external)

"""
$(TYPEDSIGNATURES)

Compute the Galerkin projection matrix of the basis `b` onto the vector Laplacian. When keyword `external=true`, 
the integral is computed over all space, assuming continuity of the poloidal field and a scalar potential in the exterior domain.
"""
diffusion(b::Basis; threads=false, external=false) = _diffusion(Val(threads), b, b; external)
diffusion(bi::Basis, bj::Basis; threads=false, external=false) = _diffusion(Val(threads), bi, bj; external)
