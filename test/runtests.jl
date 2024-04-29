using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod
using Distributed
using SparseArrays

Limace.Poly.__wiginit(100)

# include("misc.jl")
include("modes.jl")

Limace.Poly.wig_temp_free()