
"""
$(TYPEDEF)

$(TYPEDFIELDS)

This defines the Coriolis force operator for a given `basis`. The `factor` is a scalar that multiplies the operator, 
e.g. the rotation rate ``\\Omega`` or a nondimensional parameter (e.g. ``1/\\mathrm{Le}``). 
The `mat` is a sparse matrix representation of the operator, and `preassembled` indicates whether the matrix has been preassembled.

## Example usage

```julia
u = Inviscid(10)
Ω = 1.0
c = Limace.Coriolis(u, Ω)
Limace.assemble!(c)
c.mat  # sparse matrix representation of the Coriolis operator
```
"""
mutable struct Coriolis{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end


function Coriolis(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Coriolis(b, ComplexF64(factor), mat, false)
end

function assemble!(f::Coriolis; kwargs...)
    f.mat = coriolis(f.basis; kwargs...)
    f.preassembled = true
    return f.mat
end

"""
$(TYPEDSIGNATURES)

Equation (113) in [ivers_scalar_2008](@citet).
"""
@inline function C(l, m)
    return (l^2 - 1) * √((l^2 - m^2) / (4l^2 - 1))
end

##fallbacks
"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between toroidal and toroidal component, explicitly calculating the quadrature.
"""
function _coriolis_tt(bi::Ti, bj::Tj, lmna, lmnb, r, wr; Ω=2.0) where {Ti<:Basis, Tj<:Basis}
    l, m, n = lmna
    aij = _inertial_tt(bi, bj, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between poloidal and poloidal component, explicitly calculating the quadrature.
"""
function _coriolis_ss(bi::Ti, bj::Tj, lmna, lmnb, r, wr; Ω=2.0) where {Ti<:Basis, Tj<:Basis}
    l, m, n = lmna
    aij = _inertial_ss(bi, bj, lmna, lmnb, r, wr)
    return im * m * Ω / p(l) * aij
end

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between poloidal and toroidal component, explicitly calculating the quadrature. 
Following eq. (114) in [ivers_scalar_2008](@citet).
"""
function _coriolis_st(bi::Ti, bj::Tj, lmna, lmnb, r, wr; Ω=2.0) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    if lb == la-1
        _C = C(la,ma)
    elseif lb == la+1
        _C = C(la+1,ma)
    else
        return nothing
    end
    @inline _sa = r->s(Ti,bi.V,la,ma,na,r)
    @inline _tb = r->t(Tj,bj.V,lb,mb,nb,r)
    @inline f1 = r-> _∂ll(_tb,lb,la,r)
    @inline f = r -> innert(_sa,f1, la, r)
    aij = ∫dr(f,r,wr)

    return Ω / p(la) * _C * aij
end

"""
$(TYPEDSIGNATURES)

Fallback for Coriolis term between toroidal and poloidal component, explicitly calculating the quadrature.
Following eq. (112) in [ivers_scalar_2008](@citet).
"""
function _coriolis_ts(bi::Ti, bj::Tj, lmna, lmnb, r, wr; Ω=2.0) where {Ti<:Basis, Tj<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    if lb == la-1
        _C = C(la,ma)
    elseif lb == la+1
        _C = C(la+1,ma)
    else
        return nothing
    end

    @inline _ta = r->t(Ti,bi.V,la,ma,na,r)
    @inline _sb = r->s(Tj,bj.V,lb,mb,nb,r)
    @inline f1 = r-> _∂ll(_sb,lb,la,r)
    @inline f = r -> innert(_ta,f1, la, r)
    aij = ∫dr(f,r,wr)
    return Ω / p(la) * _C * aij
end

function _coriolis_poloidal_poloidal!(bi::Ti, bj::Tj, is, js, aijs, lck, lmn2k_pi, lmn2k_pj, l, m, r, wr, Ω) where {Ti<:Basis, Tj<:Basis}
    for n in nrange_p_bc(bi, l), n2 in nrange_p(bj, l)
        aij = _coriolis_ss(bi, bj, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_pi[(l, m, n)], lmn2k_pj[(l, m, n2)], aij)
    end

    return nothing
end

_coriolis_poloidal_poloidal!(b::Basis, is, js, aijs, lck, lmn2k_p, l, m, r, wr, Ω) = _coriolis_poloidal_poloidal!(b, b, is, js, aijs, lck, lmn2k_p, lmn2k_p, l, m, r, wr, Ω)

function _coriolis_poloidal_toroidal!(bi::Basis, bj::Basis, is, js, aijs, lck, _npj, lmn2k_pi, lmn2k_tj, l, l2, m, r, wr, Ω)
    for n in nrange_p_bc(bi, l), n2 in nrange_t(bj, l2)
        aij = _coriolis_st(bi, bj, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_pi[(l, m, n)], lmn2k_tj[(l2, m, n2)] + _npj, aij)
    end
    return nothing
end

_coriolis_poloidal_toroidal!(b::Basis, is, js, aijs, lck, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω) = _coriolis_poloidal_toroidal!(b, b, is, js, aijs, lck, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)

function _coriolis_toroidal_toroidal!(bi::Basis, bj::Basis, is, js, aijs, lck, _npi, _npj, lmn2k_ti, lmn2k_tj, l, m, r, wr, Ω)
    for n in nrange_t_bc(bi, l), n2 in nrange_t(bj, l)
        aij = _coriolis_tt(bi, bj, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_ti[(l, m, n)] + _npi, lmn2k_tj[(l, m, n2)] + _npj, aij)
    end
    return nothing
end

_coriolis_toroidal_toroidal!(b::Basis, is, js, aijs, lck, _np, lmn2k_t, l, m, r, wr, Ω) = _coriolis_toroidal_toroidal!(b,b, is, js, aijs, lck, _np, _np, lmn2k_t, lmn2k_t, l, m, r, wr, Ω)

function _coriolis_toroidal_poloidal!(bi::Basis, bj::Basis, is, js, aijs, lck, _npi, lmn2k_ti, lmn2k_pj, l, l2, m, r, wr, Ω)
    for n in nrange_t_bc(bi, l), n2 in nrange_p(bj, l2)
        aij = _coriolis_ts(bi, bj, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_ti[(l, m, n)] + _npi, lmn2k_pj[(l2, m, n2)], aij)
    end
    return nothing
end

_coriolis_toroidal_poloidal!(b::Basis, is, js, aijs, lck, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω) = _coriolis_toroidal_poloidal!(b,b, is, js, aijs, lck, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)

@inline function _coriolis_poloidal(bi::Basis, bj::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    lck = ReentrantLock()

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:lpmax(bi)
        for m in intersect(bi.m, -l:l)
            _coriolis_poloidal_poloidal!(bi, bj, is, js, aijs, lck, lmn2k_pi, lmn2k_pj, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l+1 > ltmax(bj)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_poloidal_toroidal!(bi, bj, is, js, aijs, lck, _npi, lmn2k_pi, lmn2k_tj, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end

_coriolis_poloidal(b::Basis; Ω::T=2.0) where T = _coriolis_poloidal(b, b; Ω)

@inline function _coriolis_toroidal(bi::Basis, bj::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    lck = ReentrantLock()

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:ltmax(bi)
        for m in intersect(bi.m, -l:l)
            _coriolis_toroidal_toroidal!(bi, bj, is, js, aijs, lck, _npi, _npj, lmn2k_ti, lmn2k_tj, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l+1 > lpmax(bj)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_toroidal_poloidal!(bi, bj, is, js, aijs, lck, _npi, lmn2k_ti, lmn2k_pj, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end

_coriolis_toroidal(b::Basis; Ω::T=2.0) where T = _coriolis_toroidal(b, b; Ω)

@inline function _coriolis_poloidal_threaded(bi::Basis, bj::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    lck = ReentrantLock()


    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    @sync for l in 1:lpmax(bi)
        for m in intersect(bi.m, -l:l)
            Threads.@spawn begin
                _coriolis_poloidal_poloidal!(bi, bj, is, js, aijs, lck, lmn2k_pi, lmn2k_pj, l, m, r, wr, Ω)
                for l2 in ((l == 1) ? (2,) : ((l+1 > ltmax(bj)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                    if l2 >= abs(m)
                        _coriolis_poloidal_toroidal!(bi, bj, is, js, aijs, lck, _npi, lmn2k_pi, lmn2k_tj, l, l2, m, r, wr, Ω)
                    end
                end
            end
        end
    end

    return vcat(is...), vcat(js...), vcat(aijs...)
end


@inline function _coriolis_toroidal_threaded(bi::Basis, bj::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lmn2k_pi = lmn2k_p_dict(bi)
    lmn2k_ti = lmn2k_t_dict(bi)
    lmn2k_pj = lmn2k_p_dict(bj)
    lmn2k_tj = lmn2k_t_dict(bj)
    _npi = np(bi)
    _npj = np(bj)
    r, wr = rquad(max(bi.N,bj.N) + 5, bi.V)
    lck = ReentrantLock()

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    @sync for l in 1:ltmax(bi)
        for m in intersect(bi.m, -l:l)
            Threads.@spawn begin
                _coriolis_toroidal_toroidal!(bi, bj, is, js, aijs, lck, _npi, _npj, lmn2k_ti, lmn2k_tj, l, m, r, wr, Ω)
                for l2 in ((l == 1) ? (2,) : ((l+1 > lpmax(bi)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                    if l2 >= abs(m)
                        _coriolis_toroidal_poloidal!(bi, bj, is, js, aijs, lck, _npi, lmn2k_ti, lmn2k_pj, l, l2, m, r, wr, Ω)
                    end
                end
            end
        end
    end

    return vcat(is...), vcat(js...), vcat(aijs...)
end


function _coriolis(::Val{false}, bi::Basis, bj::Basis; Ω::T=2.0) where T
    nui = length(bi)
    nuj = length(bj)

    is, js, aijs = _coriolis_poloidal(bi, bj; Ω)
    is2, js2, aijs2 = _coriolis_toroidal(bi, bj; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nui, nuj)
    return RHS

end

function _coriolis(::Val{true}, bi::Basis, bj::Basis; Ω::T=2.0) where T 
    nui = length(bi)
    nuj = length(bj)

    is, js, aijs = _coriolis_poloidal_threaded(bi, bj; Ω)
    is2, js2, aijs2 = _coriolis_toroidal_threaded(bi, bj; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nui, nuj)
    return RHS

end

_coriolis(::Val{true}, b::Basis; Ω::T=2.0) where T = _coriolis(Val(true), b, b; Ω)
_coriolis(::Val{false}, b::Basis; Ω::T=2.0) where T = _coriolis(Val(false), b, b; Ω)

"""
$(TYPEDSIGNATURES)

Compute the sparse Galerkin projection matrix, by projecting the basis `b` onto the Coriolis operator.
"""
coriolis(b::Basis; threads=false, Ω::T=2.0, external=false) where {T<:Number} = _coriolis(Val(threads),b, b; Ω)
coriolis(bi::Basis, bj::Basis; threads=false, Ω::T=2.0, external=false) where {T<:Number} = _coriolis(Val(threads),bi, bj; Ω)
