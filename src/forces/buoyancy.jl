"""
$(TYPEDEF)

- `ubasis::Basis{TU}`: The velocity basis for the buoyancy operator.
- `tbasis`: The temperature basis for the buoyancy operator.
- `factor::T`: : A scalar factor that multiplies the operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the operator.
- `preassembled::Bool`: A flag indicating whether the operator has been preassembled.

## Example usage

```julia
u = Inviscid(10)
T = Temperature(10)
f = Limace.Buoyancy(u,T)
Limace.assemble!(f)
f.mat = # sparse matrix representation of the buoyancy operator
```

"""
mutable struct Buoyancy{TU,TT,T} <: Forcing{2}
    ubasis::Basis{TU}
	tbasis::Basis{TT}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Buoyancy(ub::Basis, tb::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(ub), length(tb))
    return Buoyancy(ub, tb, ComplexF64(factor), mat, false)
end


function assemble!(f::Buoyancy; kwargs...)
    f.mat = buoyancy(f.ubasis, f.tbasis; kwargs...)
    f.preassembled = true
    return f.mat
end



# spherically symmetric gravity acceleration
# eq. (115) in Ivers & Phillips (2008)
function _buoyancy_sTs(::Type{TB}, ::Type{TC}, V::Volume, lmnb, lmnc, r, wr) where {TB<:Basis,TC<:Basis}
    lb,mb,nb = lmnb
    lc,mc,nc = lmnc

    @inline _Tb = r->t(TB,V,lb,mb,nb,r)
    @inline _sc = r->s(TC,V,lc,mc,nc,r)

    @inline f1 = r -> _Tb(r)

    @inline f = r-> innert(_sc,f1,lc,r)

    aij = ∫dr(f,r,wr)
    return aij
end

"""
$(TYPEDSIGNATURES)

Computes the advection term for a spherically symmetric gravity acceleration (eq. (115) in Ivers & Phillips (2008)) 
for a velocity basis `ub` and Temperature basis `tb`.
"""
function buoyancy(ub::TI, tb::TJ; kwargs...) where {TI<:Basis, TJ<:Basis}
	nu = length(ub)
	nt = length(tb)
	T = typeof(ub.V.r1)
	is,js,aijs = Int[], Int[], complex(T)[]
	lck = ReentrantLock()

    lmn2k_p_ui = lmn2k_p_dict(ub)
    lmn2k_t_tj = lmn2k_t_dict(tb)

	N = ub.N
    r,wr = rquad(N+5, ub.V)
	for li in 1:min(ltmax(tb),lpmax(ub)), mi in intersect(ub.m, -li:li)
        for ni in nrange_p_bc(ub, li)
            for nj in nrange_t(tb, li)
                aij = _buoyancy_sTs(TJ, TI, ub.V, (li, mi, nj), (li, mi, ni), r, wr)
                appendit!(is, js, aijs, lck, lmn2k_p_ui[(li,mi,ni)], lmn2k_t_tj[(li,mi,nj)], aij)
            end
        end
	end
	return sparse(is,js,aijs, nu,nt)
end