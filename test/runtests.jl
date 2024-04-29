using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod
using Distributed
using SparseArrays


function eigstarget(A, B, target; kwargs...)
    P = lu(A - target * B)
    LO = LinearMap{ComplexF64}((y, x) -> ldiv!(y, P, B * x), size(A, 2))
    pschur, history = partialschur(LO; kwargs...)
    evals, u = partialeigen(pschur)
    λ = 1 ./ evals .+ target
    return λ, u
end

Limace.Poly.__wiginit(100)

# include("misc.jl")
include("modes.jl")
include("bases_nobc.jl")

Limace.Poly.wig_temp_free()