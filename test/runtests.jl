using Limace
using Test

using LinearAlgebra
using FunctionZeros
using LinearMaps
using ArnoldiMethod
using Distributed
const DP=Limace.DiscretePart

include("pre_test.jl")
include("modes.jl")

Limace.DiscretePart.wig_temp_free()