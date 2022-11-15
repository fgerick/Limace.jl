module Limace

using LinearAlgebra
using SparseArrays
using Distributed
using DistributedArrays

# include("poly.jl")
# include("bases/viscous.jl")
# include("assemble.jl")
include("viscous.jl")
include("inviscid.jl")


end
