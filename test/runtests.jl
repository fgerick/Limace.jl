using Limace
using Test

using LinearAlgebra
using FunctionZeros
using SparseArrays

using Limace.EigenProblem: eigstarget


struct LJ22; end
Limace.Bases.s(::Type{Basis{LJ22}}, V::Volume, l, m, n, r) = r^2 * (157 - 296r^2 + 143r^4) / (16 * sqrt(182 / 3))


include("misc.jl")
include("modes.jl")
include("bases_nobc.jl")
