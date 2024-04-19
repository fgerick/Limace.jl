module Limace

using SparseArrays
using LinearAlgebra
using DocStringExtensions

include("utils.jl")
using .Utils

include("poly.jl")
using .Poly

include("quad.jl")
using .Quadrature


include("inertial.jl")
include("coriolis.jl")


include("inviscid.jl")
using .InviscidBasis

include("inviscid2.jl")
using .InviscidBasis2

include("inviscid3.jl")
using .InviscidBasis3

include("viscous.jl")
include("viscous_chen.jl")
include("insulating_allspace.jl")

include("discrete/DiscretePart.jl")

include("mhdproblem.jl")

end
