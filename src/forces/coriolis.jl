
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
Following eq. (114) in [ivers_scalar_2008](@citet).
"""
function _coriolis_st(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    if lb == la-1
        _C = C(la,ma)
    elseif lb == la+1
        _C = C(la+1,ma)
    else
        return nothing
    end
    @inline _sa = r->s(T,b.V,la,ma,na,r)
    @inline _tb = r->t(T,b.V,lb,mb,nb,r)
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
function _coriolis_ts(b::T, lmna, lmnb, r, wr; Ω=2.0) where {T<:Basis}
    la,ma,na = lmna
    lb, mb, nb = lmnb

    if lb == la-1
        _C = C(la,ma)
    elseif lb == la+1
        _C = C(la+1,ma)
    else
        return nothing
    end

    @inline _ta = r->t(T,b.V,la,ma,na,r)
    @inline _sb = r->s(T,b.V,lb,mb,nb,r)
    @inline f1 = r-> _∂ll(_sb,lb,la,r)
    @inline f = r -> innert(_ta,f1, la, r)
    aij = ∫dr(f,r,wr)
    return Ω / p(la) * _C * aij
end

function _coriolis_poloidal_poloidal!(b::T, is, js, aijs, lck, lmn2k_p, l, m, r, wr, Ω) where T<:Basis
    for n in nrange_p_bc(b, l), n2 in nrange_p(b, l)
        aij = _coriolis_ss(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_p[(l, m, n)], lmn2k_p[(l, m, n2)], aij)
    end

    return nothing
end

function _coriolis_poloidal_toroidal!(b::T, is, js, aijs, lck, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω) where T<:Basis
    for n in nrange_p_bc(b, l), n2 in nrange_t(b, l2)
        aij = _coriolis_st(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_p[(l, m, n)], lmn2k_t[(l2, m, n2)] + _np, aij)
    end
    return nothing
end


function _coriolis_toroidal_toroidal!(b::T, is, js, aijs, lck, _np, lmn2k_t, l, m, r, wr, Ω) where T<:Basis
    for n in nrange_t_bc(b, l), n2 in nrange_t(b, l)
        aij = _coriolis_tt(b, (l, m, n), (l, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_t[(l, m, n)] + _np, lmn2k_t[(l, m, n2)] + _np, aij)
    end
    return nothing
end

function _coriolis_toroidal_poloidal!(b::T, is, js, aijs, lck, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω) where T<:Basis
    for n in nrange_t_bc(b, l), n2 in nrange_p(b, l2)
        aij = _coriolis_ts(b, (l, m, n), (l2, m, n2), r, wr; Ω)
        appendit!(is, js, aijs, lck, lmn2k_t[(l, m, n)] + _np, lmn2k_p[(l2, m, n2)], aij)
    end
    return nothing
end

@inline function _coriolis_poloidal(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    lck = ReentrantLock()

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_poloidal_poloidal!(b, is, js, aijs, lck, lmn2k_p, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l+1 > ltmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_poloidal_toroidal!(b, is, js, aijs, lck, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end


@inline function _coriolis_toroidal(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], Complex{T}[]
    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)
    lck = ReentrantLock()

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            _coriolis_toroidal_toroidal!(b, is, js, aijs, lck, _np, lmn2k_t, l, m, r, wr, Ω)
            for l2 in ((l == 1) ? (2,) : ((l+1 > lpmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                if l2 >= abs(m)
                    _coriolis_toroidal_poloidal!(b, is, js, aijs, lck, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)
                end
            end
        end
    end

    return is, js, aijs
end

@inline function _coriolis_poloidal_threaded(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()


    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)


    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    @sync for l in 1:lpmax(b)
        for m in intersect(b.m, -l:l)
            Threads.@spawn begin
                _coriolis_poloidal_poloidal!(b, is, js, aijs, lck, lmn2k_p, l, m, r, wr, Ω)
                for l2 in ((l == 1) ? (2,) : ((l+1 > ltmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                    if l2 >= abs(m)
                        _coriolis_poloidal_toroidal!(b, is, js, aijs, lck, _np, lmn2k_p, lmn2k_t, l, l2, m, r, wr, Ω)
                    end
                end
            end
        end
    end

    return vcat(is...), vcat(js...), vcat(aijs...)
end


@inline function _coriolis_toroidal_threaded(b::Basis; Ω::T=2.0) where {T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_p = lmn2k_p_dict(b)
    lmn2k_t = lmn2k_t_dict(b)
    _np = np(b)
    r, wr = rquad(b.N + 5, b.V)

    #m == m2 and only l2 = l-1:l+1 needs to be considered.
    @sync for l in 1:ltmax(b)
        for m in intersect(b.m, -l:l)
            Threads.@spawn begin
                _coriolis_toroidal_toroidal!(b, is, js, aijs, lck, _np, lmn2k_t, l, m, r, wr, Ω)
                for l2 in ((l == 1) ? (2,) : ((l+1 > lpmax(b)) ? (l - 1,) : (l - 1, l + 1))) #only consider l-1 and l+1, and taking care of the upper and lower boundaries.
                    if l2 >= abs(m)
                        _coriolis_toroidal_poloidal!(b, is, js, aijs, lck, _np, lmn2k_t, lmn2k_p, l, l2, m, r, wr, Ω)
                    end
                end
            end
        end
    end

    return vcat(is...), vcat(js...), vcat(aijs...)
end


function _coriolis(::Val{false}, b::TB; Ω::T=2.0) where {TB<:Basis,T}
    nu = length(b)

    is, js, aijs = _coriolis_poloidal(b; Ω)
    is2, js2, aijs2 = _coriolis_toroidal(b; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end

function _coriolis(::Val{true}, b::TB; Ω::T=2.0) where {TB<:Basis,T}
    nu = length(b)

    is, js, aijs = _coriolis_poloidal_threaded(b; Ω)
    is2, js2, aijs2 = _coriolis_toroidal_threaded(b; Ω)

    append!(is, is2)
    append!(js, js2)
    append!(aijs, aijs2)

    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end

"""
$(TYPEDSIGNATURES)

Compute the sparse Galerkin projection matrix, by projecting the basis `b` onto the Coriolis operator.
"""
coriolis(b::Basis; threads=false, Ω::T=2.0) where {T<:Number} = _coriolis(Val(threads), b; Ω)
