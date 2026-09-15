
"""
$(TYPEDEF)

- `basis::Basis{TB}`: The basis for the inertial operator.
- `factor::T`: A scalar factor that multiplies the inertial operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the inertial operator.
- `preassembled::Bool`: A flag indicating whether the inertial operator has been preassembled.

## Example usage
```julia
u = Inviscid(10)
f = Limace.Inertial(u)
Limace.assemble!(f)
r.mat  # sparse matrix representation of the inertial operator
```
"""
mutable struct Inertial{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Inertial(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Inertial(b, ComplexF64(factor), mat, false)
end

function assemble!(f::Inertial; kwargs...)
    f.mat = sparse(inertial(f.basis; kwargs...))
    f.preassembled = true
    return f.mat
end


"""
$(TYPEDSIGNATURES)

"""
function _inertial_ss(bi::Ti, bj::Tj, lmna, lmnb, r,wr; external=false) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _sa = r->s(Ti,bi.V,la,ma,na,r)
    @inline _sb = r->s(Tj,bj.V,lb,mb,nb,r)

    @inline f = r-> inners(_sa,_sb, la, r)

    aij = ∫dr(f,r,wr)

    if external && (la==lb) && (ma==mb)
        aij += _sa(bi.V.r1)*_sb(bj.V.r1)*la^2*(la+1)
    end
    return aij
end

_inertial_ss(b::T, lmna, lmnb, r, wr; external=false) where T<:Basis = _inertial_ss(b, b, lmna, lmnb, r,wr; external)

"""
$(TYPEDSIGNATURES)

"""
function _inertial_tt(bi::Ti, bj::Tj, lmna, lmnb, r,wr) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _ta = r->t(Ti,bi.V,la,ma,na,r)
    @inline _tb = r->t(Tj,bj.V,lb,mb,nb,r)

    @inline f = r-> innert(_ta,_tb, la, r)

    aij = ∫dr(f,r,wr)
    return aij
end

_inertial_tt(b::T, lmna, lmnb, r, wr) where T<:Basis = _inertial_tt(b, b, lmna, lmnb, r,wr)


# function _inertial(::Val{false}, b::Basis; external=false)

#     T = typeof(b.V.r1)
#     is, js, aijs = Int[], Int[], complex(T)[]
#     lmn2k_p = lmn2k_p_dict(b)
#     lmn2k_t = lmn2k_t_dict(b)
#     _np = np(b)
#     r, wr = rquad(b.N + 5, b.V)
#     nu = length(b)
#     lck = ReentrantLock()

#     #m == m2 and only l==l2 needs to be considered.
#     for l in 1:lpmax(b)
#         for m in intersect(b.m, -l:l)
#             for n in nrange_p_bc(b,l), n2 in nrange_p(b,l)
#                 aij = _inertial_ss(b, (l,m,n), (l,m,n2), r,wr; external)
#                 appendit!(is, js, aijs, lck, lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
#             end
#         end
#     end

#     for l in 1:ltmax(b)
#         for m in intersect(b.m, -l:l)
#             for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
#                 aij = _inertial_tt(b, (l,m,n), (l,m,n2), r,wr)
#                 appendit!(is, js, aijs, lck, lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
#             end
#         end
#     end


#     return sparse(is, js, aijs, nu, nu)
# end

function _inertial(::Val{false}, bi::Basis, bj::Basis; external=false)

    T = promote_type(typeof(bi.V.r1),typeof(bj.V.r1))
    is, js, aijs = Int[], Int[], complex(T)[]
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    nui = length(bi)
    nuj = length(bj)
    lck = ReentrantLock()

    #m == m2 and only l==l2 needs to be considered.
    for l in 1:lpmax(bi)
        for m in intersect(bi.m, -l:l)
            for n in nrange_p_bc(bi,l), n2 in nrange_p(bj,l)
                aij = _inertial_ss(bi, bj, (l,m,n), (l,m,n2), r,wr; external)
                appendit!(is, js, aijs, lck, lmn2k_pi[(l,m,n)], lmn2k_pj[(l,m,n2)], aij)
            end
        end
    end

    for l in 1:ltmax(bi)
        for m in intersect(bi.m, -l:l)
            for n in nrange_t_bc(bi,l), n2 in nrange_t(bj,l)
                aij = _inertial_tt(bi, bj, (l,m,n), (l,m,n2), r,wr)
                appendit!(is, js, aijs, lck, lmn2k_ti[(l,m,n)] + _npi, lmn2k_tj[(l,m,n2)] + _npj, aij)
            end
        end
    end


    return sparse(is, js, aijs, nui, nuj)
end

_inertial(::Val{false}, b::Basis; external=false) = _inertial(Val(false), b, b; external)

function _inertial(::Val{true}, bi::Basis, bj::Basis; external=false)

    T = promote_type(typeof(bi.V.r1),typeof(bj.V.r1))
    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    nui = length(bi)
    nuj = length(bj)

    #m == m2 and only l==l2 needs to be considered.
    @sync begin
        for l in 1:lpmax(bi)
            for m in intersect(bi.m, -l:l)
                for n in nrange_p_bc(bi,l), n2 in nrange_p(bj,l)
                    Threads.@spawn begin
                        aij = _inertial_ss(bi,bj, (l,m,n), (l,m,n2), r,wr; external)
                        appendit!(is, js, aijs, lck, lmn2k_pi[(l,m,n)], lmn2k_pj[(l,m,n2)], aij)
                    end
                end
            end
        end

        for l in 1:ltmax(bi)
            for m in intersect(bi.m, -l:l)
                for n in nrange_t_bc(bi,l), n2 in nrange_t(bj,l)
                    Threads.@spawn begin
                        aij = _inertial_tt(bi, bj, (l,m,n), (l,m,n2), r,wr)
                        appendit!(is, js, aijs, lck, lmn2k_ti[(l,m,n)] + _npi, lmn2k_tj[(l,m,n2)] + _npj, aij)
                    end
                end
            end
        end
    end


    return  sparse(vcat(is...), vcat(js...), vcat(aijs...), nui, nuj)
end

_inertial(::Val{true}, b::Basis; external=false) = _inertial(Val(true), b, b; external)

"""
$(TYPEDSIGNATURES)

Compute the Galerkin projection matrix of basis `b` onto itself, i.e. the inner products. 
Also known as the mass matrix.
"""
inertial(b::Basis; threads=false, external=false) = inertial(b, b; threads, external)

inertial(bi::Basis, bj::Basis; threads=false, external=false) = _inertial(Val(threads), bi, bj; external)
