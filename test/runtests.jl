using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod
using Distributed
const DP=Limace.DiscretePart

Limace.DiscretePart.__wiginit(100)

# include("misc.jl")
# include("pre_test.jl")
include("modes.jl")

Limace.DiscretePart.wig_temp_free()