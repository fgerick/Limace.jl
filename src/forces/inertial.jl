
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


function _inertial(::Val{false}, b::Basis; external=false)

    is, js, aijs = Int[], Int[], Complex{Float64}[]
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

function _inertial(::Val{true}, b::Basis; external=false)

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [Complex{Float64}[] for _ in 1:_nt]

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    nu = length(b)

    #m == m2 and only l==l2 needs to be considered.
    @sync begin
        for l in 1:lpmax(b)
            for m in intersect(b.m, -l:l)
                for n in nrange_p_bc(b,l), n2 in nrange_p(b,l)
                    Threads.@spawn begin
                        id = Threads.threadid()
                        aij = _inertial_ss(b, (l,m,n), (l,m,n2), r,wr; external)
                        appendit!(is[id], js[id], aijs[id], lmn2k_p[(l,m,n)], lmn2k_p[(l,m,n2)], aij)
                    end
                end
            end
        end

        for l in 1:ltmax(b)
            for m in intersect(b.m, -l:l)
                for n in nrange_t_bc(b,l), n2 in nrange_t(b,l)
                    Threads.@spawn begin
                        id = Threads.threadid()
                        aij = _inertial_tt(b, (l,m,n), (l,m,n2), r,wr)
                        appendit!(is[id], js[id], aijs[id], lmn2k_t[(l,m,n)] + _np, lmn2k_t[(l,m,n2)] + _np, aij)
                    end
                end
            end
        end
    end


    return  sparse(vcat(is...), vcat(js...), vcat(aijs...), nu, nu)
end

"""
$(TYPEDSIGNATURES)

Compute the Galerkin projection matrix of basis `b` onto itself, i.e. the inner products. 
Also known as the mass matrix.
"""
inertial(b::Basis; threads=false, external=false) = _inertial(Val(threads), b; external)
