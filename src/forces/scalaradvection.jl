# 𝐮⋅∇T
# Equation (118) in Ivers & Phillips (2008)

"""
$(TYPEDEF)

- `tbasis::Basis{TT}`: The temperature basis for the scalar advection operator.
- `ubasis::Basis{TU}`: The basis for the scalar advection operator.
- `T0`: The background temperature, which can be a single `BasisElement` or a collection of them.
- `factor::T`: : A scalar factor that multiplies the operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the operator.
- `preassembled::Bool`: A flag indicating whether the operator has been preassembled.

## Example usage

```julia
u = Limace.Inviscid(10)
T = Limace.Temperature(10)
T0 = BasisElement(t, Toroidal, (1,0,1))
f = Limace.ScalarAdvectionT0(T,u,T0)
Limace.assemble!(f)
f.mat # sparse matrix representation of the scalar advection operator
```

"""
mutable struct ScalarAdvectionT0{TT,TU,T} <: Forcing{2}
    tbasis::Basis{TT}
    ubasis::Basis{TU}
    T0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function ScalarAdvectionT0(tb::Basis, ub::Basis, T0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(tb), length(ub))
    return ScalarAdvectionT0(tb, ub, T0, ComplexF64(factor), mat, false)
end


function assemble!(f::ScalarAdvectionT0; kwargs...)
    f.mat = sum(scalaradvection(f.tbasis, f.ubasis, T0; kwargs...) for T0 in f.T0)
    f.preassembled = true
    return f.mat
end


"""
$(TYPEDEF)

- `tbasis::Basis{TT}`: The temperature basis for the scalar advection operator.
- `U0`: The background flow field, which can be a single `BasisElement` or a collection of them.
- `factor::T`: : A scalar factor that multiplies the operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the operator.
- `preassembled::Bool`: A flag indicating whether the operator has been preassembled.

## Example usage

```julia
u = Limace.Inviscid(10)
T = Limace.Temperature(10)
U0 = BasisElement(u, Toroidal, (1,0,1))
f = Limace.ScalarAdvectionU0(T,U0)
Limace.assemble!(f)
f.mat # sparse matrix representation of the scalar advection operator
```

"""
mutable struct ScalarAdvectionU0{TT,T} <: Forcing{1}
    tbasis::Basis{TT}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function ScalarAdvectionU0(tb::Basis, U0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(tb), length(tb))
    return ScalarAdvectionU0(tb, U0, ComplexF64(factor), mat, false)
end


function assemble!(f::ScalarAdvectionU0; kwargs...)
    f.mat = sum(scalaradvection(f.tbasis, U0, f.tbasis; kwargs...) for U0 in f.U0)
    f.preassembled = true
    return f.mat
end

# Aabc
function _scalaradvection_sTT(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _sa = r->s(TA,V,la,ma,na,r)
    @inline _Tb = r->t(TB,V,lb,mb,nb,r)
    @inline _Tc = r->t(TC,V,lc,mc,nc,r)

    @inline f1 = r -> (r*p(la)*_sa(r)*∂(_Tb,r) + 1/2*(p(la)+p(lb)-p(lc))*∂(r->r*_sa(r),r)*_Tb(r))/(p(lc)*r^2)

    @inline f = r-> -innert(_Tc,f1,lc,r)

    aij = ∫dr(f,r,wr)
    return aij
end

# Eabc 
function _scalaradvection_tTT(::Type{TA}, ::Type{TB}, ::Type{TC}, V::Volume, lmna, lmnb, lmnc, r, wr) where {TA<:Basis,TB<:Basis,TC<:Basis}
    la,ma,na = lmna
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _ta = r->t(TA,V,la,ma,na,r)
    @inline _Tb = r->t(TB,V,lb,mb,nb,r)
    @inline _Tc = r->t(TC,V,lc,mc,nc,r)

    @inline f1 = r -> -_ta(r)*_Tb(r)/(p(lc)*r)

    @inline f = r-> -innert(_Tc,f1,lc,r)

    aij = ∫dr(f,r,wr)
    return aij
end

function _scalaradvection(::Val{false}, bti::TI, buj::TJ, t0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_t_ti = lmn2k_t_dict(bti)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = t0.lmn
    @assert bti.N == buj.N "Use same resolution for bases!"
    N = bti.N
    rwrs = [rquad(n + l0 + n0 + 5, bti.V) for n in 1:N]

    npu = length(lmn2k_p_uj)

    for li in 1:ltmax(bti), mi in intersect(bti.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
            A = adamgaunt(lj,l0,li, mj, m0, mi)
            _crossterm!(bti, buj, t0, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_p_uj, nrange_t_bc, nrange_p, _scalaradvection_sTT, A)
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            E = elsasser(lj, l0, li, mj, m0, mi)
            _crossterm!(bti, buj, t0, is, js, aijs, lck, 0, npu, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_t_uj, nrange_t_bc, nrange_t, _scalaradvection_tTT, E)
        end
    end

    nmatt = length(bti)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatt, nmatu)
end

function _scalaradvection(::Val{true}, bti::TI, buj::TJ, t0::BasisElement{T0,Toroidal,T}) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

    is, js, aijs = Int[], Int[], complex(T)[]
    lck = ReentrantLock()

    lmn2k_t_ti = lmn2k_t_dict(bti)

    lmn2k_p_uj = lmn2k_p_dict(buj)
    lmn2k_t_uj = lmn2k_t_dict(buj)

    l0, m0, n0 = t0.lmn
    @assert bti.N == buj.N "Use same resolution for bases!"
    N = bti.N
    rwrs = [rquad(n + l0 + n0 + 5, bti.V) for n in 1:N]

    npu = length(lmn2k_p_uj)

    @sync for li in 1:ltmax(bti), mi in intersect(bti.m, -li:li)
        mj = adamgaunt_mjs(mi, m0)
        for lj in adamgaunt_ljs(li, l0, mj, lpmax(buj))
             Threads.@spawn begin
                A = adamgaunt(lj,l0,li, mj, m0, mi)
                _crossterm!(bti, buj, t0, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_p_uj, nrange_t_bc, nrange_p, _scalaradvection_sTT, A)
             end
        end
        mj = elsasser_mjs(mi, m0)
        for lj in elsasser_ljs(li, l0, mj, ltmax(buj))
            Threads.@spawn begin
                E = elsasser(lj, l0, li, mj, m0, mi)
                _crossterm!(bti, buj, t0, is, js, aijs, lck, 0, npu, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_t_uj, nrange_t_bc, nrange_t, _scalaradvection_tTT, E)
            end
        end
    end

    nmatt = length(bti)
    nmatu = length(buj)

    return sparse(is, js, aijs, nmatt, nmatu)
end

# function _scalaradvection(::Val{false}, bti::TI, U0::BasisElement{T0,Toroidal,T}, btj::TJ) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

#     is, js, aijs = Int[], Int[], complex(T)[]
#     lck = ReentrantLock()

#     lmn2k_t_ti = lmn2k_t_dict(bti)
#     lmn2k_t_tj = lmn2k_t_dict(btj)

#     l0, m0, n0 = U0.lmn
#     @assert bti.N == btj.N "Use same resolution for bases!"
#     N = bti.N
#     rwrs = [rquad(n + l0 + n0 + 1, bti.V) for n in 1:N]

#     for li in 1:ltmax(bti), mi in intersect(bti.m, -li:li)
#         mj = elsasser_mjs(mi, m0)
#         for lj in elsasser_ljs(li, l0, mj, ltmax(btj))
#             E = elsasser(l0, lj, li, m0, mj, mi)
#             _crossterm!(bti, U0, btj, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_t_tj, nrange_t_bc, nrange_t, _scalaradvection_tTT, E)
#         end
#     end

#     nmatti = length(bti)
#     nmattj = length(btj)

#     return sparse(is, js, aijs, nmatti, nmattj)
# end

# function _scalaradvection(::Val{false}, bti::TI, U0::BasisElement{T0,Poloidal,T}, btj::TJ) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}

#     is, js, aijs = Int[], Int[], complex(T)[]
#     lck = ReentrantLock()

#     lmn2k_t_ti = lmn2k_t_dict(bti)
#     lmn2k_t_tj = lmn2k_t_dict(btj)

#     l0, m0, n0 = U0.lmn
#     @assert bti.N == btj.N "Use same resolution for bases!"
#     N = bti.N
#     rwrs = [rquad(n + l0 + n0 + 1, bti.V) for n in 1:N]

#     for li in 1:ltmax(bti), mi in intersect(bti.m, -li:li)
# 		mj = adamgaunt_mjs(mi, m0)
#         for lj in adamgaunt_ljs(li, l0, mj, ltmax(btj))
#             A = adamgaunt(l0,lj,li, m0, mj, mi)
#             _crossterm!(bti, U0, btj, is, js, aijs, lck, 0, 0, li, mi, lj, mj, rwrs, lmn2k_t_ti, lmn2k_t_tj, nrange_t_bc, nrange_t, _scalaradvection_sTT, A)
#         end
#     end

#     nmatti = length(bti)
#     nmattj = length(btj)

#     return sparse(is, js, aijs, nmatti, nmattj)
# end


"""
$(TYPEDSIGNATURES)

Computes the scalar advection term for a background temperature `t0`, a temperature basis `bti` and a velocity basis `buj`.
Currently, scalar fields are implemented as being `Toroidal` scalars (may change in the future).
"""
function scalaradvection(bti::TI, buj::TJ, t0::BasisElement{T0,Toroidal,T}; threads=false, external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,T}
    return _scalaradvection(Val(threads), bti, buj, t0)
end

# """
# $(TYPEDSIGNATURES)

# Computes the scalar advection term for a poloidal/toroidal background flow `U0`, a temperature basis `bti` and a temperature basis `btj`.
# """
# function scalaradvection(bti::TI, U0::BasisElement{T0,TH,T}, btj::TJ; threads=false, external=true) where {TI<:Basis,TJ<:Basis,T0<:Basis,TH<:Helmholtz,T}
#     return _scalaradvection(Val(threads), bti, U0, btj)
# end



# function _induction_new(::Val{false}, bti::Basis{TI,V}, buj::Basis{TJ,V}, B0::BasisElement{Basis{T0,V},Poloidal,T}; external=true) where {TI,TJ,T0,T, V<:Volume}

#     is, js, aijs = Int[], Int[], complex(T)[]
#     lck = ReentrantLock()

#     lmn2k_p_bi = lmn2k_p_dict(bti)
#     lmn2k_t_bi = lmn2k_t_dict(bti)

#     lmn2k_p_uj = lmn2k_p_dict(buj)
#     lmn2k_t_uj = lmn2k_t_dict(buj)

#     l0, m0, n0 = B0.lmn
#     @assert bti.N == buj.N "Use same resolution for bases!"
#     N = bti.N
#     rwrs = [rquad(n + l0 + n0 + 1, bti.V) for n in 1:N]

#     npb = length(lmn2k_p_bi)
#     npu = length(lmn2k_p_uj)
    
#     for li in 1:lpmax(bti)
#         for ni in nrange_p_bc(bti, li)
#             for lj in adamgaunt_ljs(li, l0, 0, lpmax(buj))
#                _crossterm_m_adamgaunt!(bti,buj,B0, is, js, aijs, lck, 0,0, li,ni,lj, rwrs, lmn2k_p_bi, lmn2k_p_uj, nrange_p, lpmax, _induction_sSS; external)  
#             end
#             for lj in elsasser_ljs(li, l0, 0, ltmax(buj))
#                _crossterm_m_elsasser!(bti,buj,B0, is, js, aijs, lck, 0,npu, li,ni,lj, rwrs, lmn2k_p_bi, lmn2k_t_uj, nrange_t, ltmax, _induction_tSS; external)  
#             end
#         end
#     end

#     for li in 1:ltmax(bti)
#         for ni in nrange_t_bc(bti, li)
#             for lj in  elsasser_ljs(li, l0, 0, lpmax(buj))
#                _crossterm_m_elsasser!(bti,buj,B0, is, js, aijs, lck, npb,0, li,ni,lj, rwrs, lmn2k_t_bi, lmn2k_p_uj, nrange_p, lpmax, _induction_sST)  
#             end
#             for lj in  adamgaunt_ljs(li, l0, 0, ltmax(buj))
#                _crossterm_m_adamgaunt!(bti, buj, B0, is, js, aijs, lck, npb, npu, li,ni,lj, rwrs, lmn2k_t_bi, lmn2k_t_uj, nrange_t, ltmax, _induction_tST)
#             end
#         end
#     end


#     nmatb = length(bti)
#     nmatu = length(buj)

#     return sparse(is, js, aijs, nmatb, nmatu)
# end
