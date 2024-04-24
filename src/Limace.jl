module Limace

using SparseArrays
using LinearAlgebra
using DocStringExtensions

include("Utils.jl")
using .Utils

include("poly.jl")
using .Poly

include("quad.jl")
using .Quadrature

include("Bases.jl")
using .Bases

# forces
include("forces/inertial.jl")
include("forces/coriolis.jl")
include("forces/diffusion.jl")


include("bases/Inviscid.jl")
using .InviscidBasis

include("bases/inviscid2.jl")
using .InviscidBasis2

include("bases/InviscidNoBC.jl")
using .InviscidBasisNoBC

include("bases/Insulating.jl")
using .InsulatingBasis

# include("bases/viscous.jl")
# include("bases/viscous_chen.jl")

# include("discrete/DiscretePart.jl")

# include("mhdproblem.jl")

end
