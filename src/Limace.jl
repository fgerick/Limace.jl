module Limace

using LinearAlgebra
using SparseArrays
using SpecialFunctions
# using Distributed
# using DistributedArrays
using Random

# include("poly.jl")
# include("bases/viscous.jl")
# include("assemble.jl")
include("viscous.jl")
include("inviscid.jl")
include("viscous_chen.jl")
include("insulating_allspace.jl")

include("discrete/DiscretePart.jl")

include("mhdproblem.jl")

end
