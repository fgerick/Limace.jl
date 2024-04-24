module DiscretePart

# using Limace
using ForwardDiff
using SparseArrays
using SpecialFunctions
using Wigxjpf
using FastGaussQuadrature
using Distributed
using DistributedArrays
using Random
using DocStringExtensions

using ..Utils
import ..InviscidBasis
import ..InsulatingMFBasis
using ..Poly


# include("quad.jl")
include("pt_scalars.jl")
include("pt_scalars_pre.jl")

include("induction.jl")
include("lorentz.jl")
include("lorentz_pre.jl")
include("advection_pre.jl")
include("induction_pre.jl")
include("diffusion.jl")
# include("inertial.jl")
# include("coriolis.jl")


end