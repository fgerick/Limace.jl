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

#only for dev
function __wiginit(N)
    wig_table_init(2N, 9)
    wig_temp_init(2N)
end

function __wiginit_thread(N)
    wig_table_init(2N, 9)
    wig_thread_temp_init(2N)
end

include("quad.jl")
include("pt_scalars.jl")
include("pt_scalars_pre.jl")

include("induction.jl")
include("lorentz.jl")
include("lorentz_pre.jl")
include("advection_pre.jl")
include("induction_pre.jl")
include("diffusion.jl")
include("inertial.jl")
include("coriolis.jl")


end