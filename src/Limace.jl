module Limace


include("utils.jl")
include("poly.jl")

include("viscous.jl")
include("inviscid.jl")
include("viscous_chen.jl")
include("insulating_allspace.jl")

include("discrete/DiscretePart.jl")

include("mhdproblem.jl")

using .Utils, .Poly
using .InviscidBasis

end
