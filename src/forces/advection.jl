"""
$(TYPEDEF)

- `basis::Basis{TU}`: The basis for the advection operator.
- `U0`: The background flow, which can be a single `BasisElement` or a collection of them.
- `factor::T`: A scalar factor that multiplies the advection operator, defaulting to `1.0`.
- `mat::SparseMatrixCSC{ComplexF64}`: A sparse matrix representation of the advection operator.
- `preassembled::Bool`: A flag indicating whether the advection operator has been preassembled.

## Example usage
```julia
u = Inviscid(10)
U0 = BasisElement(u, Toroidal, (1,0,0))
a = Limace.Advection(u, U0)
Limace.assemble!(a)
a.mat  # sparse matrix representation of the advection operator
```
"""
mutable struct Advection{TU,T} <: Forcing{1}
    basis::Basis{TU}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Advection(ub::Basis, U0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(ub), length(ub))
    return Advection(ub, U0, ComplexF64(factor), mat, false)
end


function assemble!(f::Advection; kwargs...)
    f.mat = sum(advection(f.basis, U0; kwargs...) for U0 in f.U0)
    f.preassembled = true
    return f.mat
end

"""
$(TYPEDSIGNATURES)

Computes the advection term for a poloidal/toroidal background flow `U0`, a velocity basis `u`. Equivalent to `-lorentz(u,u,U0)`.
"""
advection(u::Tu, U0::BasisElement{T0,TH,T}; threads=false) where {Tu<:Basis,T0<:Basis,TH<:Helmholtz,T} = -_lorentz(Val(threads), u, u, U0)
