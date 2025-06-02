"""
$(TYPEDEF)

- `basis::Basis{TB}`: The basis for the induction operator.
- `U0`: The background flow, which can be a single `BasisElement` or a collection of them.
- `factor::T`: : A scalar factor that multiplies the induction operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the induction operator.
- `preassembled::Bool`: A flag indicating whether the induction operator has been preassembled.

## Example usage

```julia
u = Inviscid(10)
b = Insulating(10)
U0 = BasisElement(u, Toroidal, (2,0,1))
f = Limace.InductionU0(b, U0)
Limace.assemble!(f)
f.mat # sparse matrix representation of the induction operator
```
"""
mutable struct InductionU0{TB,T} <: Forcing{1}
    basis::Basis{TB}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

- `bbasis::Basis{TB}`: The magnetic field basis for the induction operator.
- `ubasis::Basis{TB}`: The velocity basis for the induction operator.
- `B0`: The background magnetic field, which can be a single `BasisElement` or a collection of them.
- `factor::T`: : A scalar factor that multiplies the induction operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the induction operator.
- `preassembled::Bool`: A flag indicating whether the induction operator has been preassembled.

```julia
u = Inviscid(10)
b = Insulating(10)
B0 = BasisElement(b, Poloidal, (2,0,1))
f = Limace.InductionB0(b, u, B0)
Limace.assemble!(f)
f.mat # sparse matrix representation of the induction operator
```
"""
mutable struct InductionB0{TB,TU,T} <: Forcing{2}
    bbasis::Basis{TB}
    ubasis::Basis{TU}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function InductionU0(b::Basis, U0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return InductionU0(b, U0, ComplexF64(factor), mat, false)
end

function InductionB0(bb::Basis, ub::Basis, B0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(bb), length(ub))
    return InductionB0(bb, ub, B0, ComplexF64(factor), mat, false)
end


function assemble!(f::InductionU0; kwargs...)
    f.mat = sum(induction(f.basis, U0, f.basis; kwargs...) for U0 in f.U0)
    f.preassembled = true
    return f.mat
end

function assemble!(f::InductionB0; kwargs...)
    f.mat = sum(induction(f.bbasis, f.ubasis, B0; kwargs...) for B0 in f.B0)
    f.preassembled = true
    return f.mat
end

##
## poloidal induction equation
##

"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(pⱼ×Sₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_sSS(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> (-p(la) * (-p(la) + p(lb) + p(lc)) * s(TA, V, la, ma, na, r) * ∂(r -> r * s(TB, V, lb, mb, nb, r), r) +
                       p(lb) * (p(la) - p(lb) + p(lc)) * s(TB, V, lb, mb, nb, r) * ∂(r -> r * s(TA, V, la, ma, na, r), r)) / (2r^2 * p(lc))

    @inline f = r -> inners(f1, r -> s(TC, V, lc, mc, nc, r), lc, r)

    aij = ∫dr(f, r, wr) 

    if external
        aij += f1(V.r1) * s(TC, V, lc, mc, nc, V.r1) * p(lc) * lc 
    end
    return aij
end


#poloidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(pⱼ×Tₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_sTS(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> p(la) * s(TA, V, la, ma, na, r) * t(TB, V, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, V, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end


#toroidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Sᵢ⋅∇×(qⱼ×Sₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_tSS(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr; external=true) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline f1 = r -> p(lb) * t(TA, V, la, ma, na, r) * s(TB, V, lb, mb, nb, r) / (r * p(lc))
    @inline f = r -> inners(r -> s(TC, V, lc, mc, nc, r), f1, lc, r)

    aij = ∫dr(f, r, wr) 

    #add contribution from external ∫dV (1<r<∞), 
    #if toroidal velocity is not 0 at r=11 

    if external
        aij += f1(V.r1) * s(TC,V, lc, mc, nc, V.r1) * lc * p(lc) 
    end

    return aij
end


#toroidal flow, toroidal B0 
#always 0


##
## toroidal induction equation
##

#poloidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(pⱼ×Sₖ) dV. Surface integral is given by Elsasser variable.
"""
function _induction_sST(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _sa = r -> s(TA, V, la, ma, na, r)
    @inline _Sb = r -> s(TB, V, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, V, lc, mc, nc, r)

    @inline f1 = r -> ((p(la) + p(lb) + p(lc)) * _sa(r) * _Sb(r) -
                       (p(la) + p(lb) - p(lc)) * (r * ∂(_sa, r) * _Sb(r) + r * _sa(r) * ∂(_Sb, r) + r^2 * ∂(_sa, r) * ∂(_Sb, r)) -
                       p(la) * r^2 * _sa(r) * ∂(r -> ∂(_Sb, r), r) - p(lb) * r^2 * _Sb(r) * ∂(r -> ∂(_sa, r), r)) / (r^3 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end


#poloidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(pⱼ×Tₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_sTT(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _sa = r -> s(TA, V, la, ma, na, r)
    @inline _Tb = r -> t(TB, V, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, V, lc, mc, nc, r)

    @inline f1 = r -> (-p(lc) * (p(la) + p(lb) - p(lc)) * (_sa(r) * _Tb(r) + r * ∂(_sa, r) * _Tb(r)) +
                       p(la) * (p(la) - p(lb) - p(lc)) * (r * ∂(_sa, r) * _Tb(r) + r * _sa(r) * ∂(_Tb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end

#toroidal flow, poloidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(qⱼ×Sₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_tST(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _ta = r -> t(TA, V, la, ma, na, r)
    @inline _Sb = r -> s(TB, V, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, V, lc, mc, nc, r)

    @inline f1 = r -> (p(lc) * (p(la) + p(lb) - p(lc)) * (_ta(r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r)) -
                       p(lb) * (p(lb) - p(la) - p(lc)) * (r * ∂(_ta, r) * _Sb(r) + r * _ta(r) * ∂(_Sb, r))) / (2r^2 * p(lc))

    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end

#toroidal flow, toroidal B0
"""
$(TYPEDSIGNATURES)

Computes radial integral of ∫Tᵢ⋅∇×(qⱼ×Tₖ) dV. Surface integral is given by Adam-Gaunt variable.
"""
function _induction_tTT(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la, ma, na = lmna
    lb, mb, nb = lmnb
    lc, mc, nc = lmnc

    @inline _ta = r -> t(TA, V, la, ma, na, r)
    @inline _Tb = r -> t(TB, V, lb, mb, nb, r)
    @inline _Tc = r -> t(TC, V, lc, mc, nc, r)

    @inline f1 = r -> _ta(r) * _Tb(r) / r
    @inline f = r -> innert(_Tc, f1, lc, r)

    aij = ∫dr(f, r, wr) 
    return aij
end



#matrix assembly

"""
$(TYPEDSIGNATURES)

For a combination of `mi`, `m0` and `mj` only `mj=mi-m0` is nonzero for the Adam-Gaunt variables.
"""
@inline adamgaunt_mjs(mi, m0) = mi - m0

"""
$(TYPEDSIGNATURES)

For a combination of `li`, `l0` and `lj` only the even `lj` that satisfy `|l0-li|≤lj≤l0+li` are nonzero for the Adam-Gaunt variables.
"""
@inline function adamgaunt_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if iseven(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end

"""
$(TYPEDSIGNATURES)

For a combination of `mi`, `m0` and `mj` only `mj=mi-m0` is nonzero for the Elsasser variables.
"""
@inline elsasser_mjs(mi, m0) = mi - m0

"""
$(TYPEDSIGNATURES)

For a combination of `li`, `l0` and `lj` only the odd `lj` that satisfy `|l0-li|≤lj≤l0+li` are nonzero for the Elsasser variables.
"""
@inline function elsasser_ljs(li, l0, mj, lmax)
    lj0 = max(abs(l0 - li), max(1, abs(mj)))
    ljmax = min(li+l0, lmax)
    if isodd(li + l0 + lj0)
        return lj0:2:ljmax
    else
        return lj0+1:2:ljmax
    end
end

"""
$(TYPEDSIGNATURES)

Fallback functions for `_crossterm` term for `U0`. Write specialized function to include e.g. bandedness in `n`.
"""
@inline function _crossterm!(bi::TI, B0::BasisElement{T0,PT,T}, bj::TJ, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(max(bi.N,bj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(T0, TJ, TI, bi.V, B0.lmn, lmnj, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Fallback functions for `_crossterm!` term for `B0`. Write specialized function to include e.g. bandedness in `n`.
"""
@inline function _crossterm!(bi::TI, bj::TJ, B0::BasisElement{T0,PT,T}, is, js, aijs, i0, j0,
    li, mi, lj, mj, rwrs, lmn2k_bi, lmn2k_bj, nrangefi, nrangefj, indf, EA; kwargs...) where {TI<:Basis,TJ<:Basis,T0<:Basis,PT<:Helmholtz,T}
    l0,m0,n0 = B0.lmn
    for ni in nrangefi(bi, li)
        for nj in nrangefj(bj, lj)
            r, wr = rwrs[min(max(bi.N,bj.N), li ÷ 2 + ni + lj ÷ 2 + nj + 1 + l0 + n0)]
            lmni = (li, mi, ni)
            lmnj = (lj, mj, nj)
            aij = indf(TJ, T0, TI, bi.V, lmnj, B0.lmn, lmni, r, wr; kwargs...)*EA
            appendit!(is, js, aijs, lmn2k_bi[lmni] + i0, lmn2k_bj[lmnj] + j0, aij*B0.factor)
        end
    end
    return nothing
end



function _induction(::Val{false}, bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sSS,A; external)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj,nrange_p_bc, nrange_t, _induction_tSS, E; external)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj,nrange_t_bc, nrange_t, _induction_tST, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj,nrange_t_bc, nrange_p, _induction_sST, E)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end


function _induction(::Val{false}, bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 5, bbi.V) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sTS, E)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sTT, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bbi, buj, B0, is, js, aijs, npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tTT, E)
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatb, nmatu)
end


function _induction(::Val{false}, bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_sSS, A; external)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, 0, npbj, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, nrange_p_bc, nrange_t, _induction_sTS, E)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_sTT, A)
        end
        mj = elsasser_mjs(mi,m0)
        for lj in elsasser_ljs(li,l0,mj, lpmax(bbj))
            E = elsasser(l0,lj,li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_sST, E)
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(is, js, aijs, nmatbi, nmatbj)
end


function _induction(::Val{false}, bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_tSS, E; external)
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_tST, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            _crossterm!(bbi, U0, bbj, is, js, aijs, npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_tTT, E)
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(is, js, aijs, nmatbi, nmatbj)
end

function _induction(::Val{true}, bbi::TI, buj::TJ, B0::BasisElement{T0,Poloidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]

    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)

    @sync begin
    for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                id = Threads.threadid()
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sSS, A; external)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                id = Threads.threadid()
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], 0, npu, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_uj, nrange_p_bc, nrange_t, _induction_tSS, E; external)
            end
        end
    end

    for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(buj))
            Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                id = Threads.threadid()
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tST, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                id = Threads.threadid()
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sST, E)
            end
        end
    end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
end

function _induction(::Val{true}, bbi::TI, buj::TJ, B0::BasisElement{T0,Toroidal,T}; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = B0.lmn
    @assert bbi.N == buj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npb = length(lmn2k_p_bi)
    npu = length(lmn2k_p_uj)
    @sync for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(buj))
            Threads.@spawn begin
                id = Threads.threadid()
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p_bc, nrange_p, _induction_sTS, E)
            end
        end
    end

    @sync for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            Threads.@spawn begin
                id = Threads.threadid()
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], npb, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_t_bc, nrange_p, _induction_sTT, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            Threads.@spawn begin
                id = Threads.threadid()
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bbi, buj, B0, is[id], js[id], aijs[id], npb, npu, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t_bc, nrange_t, _induction_tTT, E)
            end
        end
    end

    nmatb = length(bbi)
    nmatu = length(buj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatb, nmatu)
end

function _induction(::Val{true}, bbi::TI, U0::BasisElement{T0,Poloidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    @sync for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_sSS, A; external)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, npbj, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_t_bj, nrange_p_bc, nrange_t, _induction_sTS, E)
            end
        end
    end

    @sync for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, ltmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_sTT, A)
            end
        end
        mj = elsasser_mjs(mi,m0)
        for lj in elsasser_ljs(li,l0,mj, lpmax(bbj))
            E = elsasser(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_sST, E)
            end
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatbi, nmatbj)
end

function _induction(::Val{true}, bbi::TI, U0::BasisElement{T0,Toroidal,T}, bbj::TJ; external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    
    _nt = Threads.nthreads()
    is, js, aijs = [Int[] for _ in 1:_nt], [Int[] for _ in 1:_nt], [complex(T)[] for _ in 1:_nt]


    lmn2k_p_bi = lmn2k_p_dict(bbi)
    lmn2k_t_bi = lmn2k_t_dict(bbi)

    lmn2k_p_bj = lmn2k_p_dict(bbj)
    lmn2k_t_bj = lmn2k_t_dict(bbj)

    l0, m0, n0 = U0.lmn
    @assert bbi.N == bbj.N "Use same resolution for bases!"
    N = bbi.N
    rwrs = [rquad(n + l0 + n0 + 1, bbi.V) for n in 1:N]

    npbi = length(lmn2k_p_bi)
    npbj = length(lmn2k_p_bj)
    @sync for li in 1:lpmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, lpmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], 0, 0, li, mi, lj, mj, rwrs, lmn2k_p_bi, lmn2k_p_bj, nrange_p_bc, nrange_p, _induction_tSS, E; external)
            end
        end
    end

    @sync for li in 1:ltmax(bbi), mi in intersect(bbi.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(bbj))
            A = adamgaunt(l0,lj,li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, 0, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_p_bj, nrange_t_bc, nrange_p, _induction_tST, A)
            end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(bbj))
            E = elsasser(l0, lj, li, m0, mj, mi)
            Threads.@spawn begin
                id = Threads.threadid()
                _crossterm!(bbi, U0, bbj, is[id], js[id], aijs[id], npbi, npbj, li, mi, lj, mj, rwrs, lmn2k_t_bi, lmn2k_t_bj, nrange_t_bc, nrange_t, _induction_tTT, E)
            end
        end
    end

    nmatbi = length(bbi)
    nmatbj = length(bbj)

    return sparse(vcat(is...), vcat(js...), vcat(aijs...), nmatbi, nmatbj)
end

"""
$(TYPEDSIGNATURES)

Computes the induction term for a poloidal/toroidal background magnetic field `B0`, a magnetic field basis `bbi` and a velocity basis `buj`.
"""
induction(bbi::TI, buj::TJ, B0::BasisElement{T0,TH,T}; threads=false, external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,TH<:Helmholtz,T} = _induction(Val(threads), bbi, buj, B0; external)

"""
$(TYPEDSIGNATURES)

Computes the induction term for a poloidal/toroidal background velocity `U0`, a magnetic field basis `bbi` and a magnetic field basis `bbj`.
"""
induction(bbi::TI, U0::BasisElement{T0,TH,T}, bbj::TJ; threads=false, external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,TH<:Helmholtz,T} = _induction(Val(threads), bbi, U0, bbj; external)
