module DiscretePart

using Limace
using ForwardDiff
using SparseArrays
using SpecialFunctions
using Wigxjpf
using FastGaussQuadrature
using Distributed
using DistributedArrays
using Random

#only for dev
function __wiginit(N)
    wig_table_init(2N, 9)
    wig_temp_init(2N)
end

include("poly.jl")
include("quad.jl")
include("pt_scalars.jl")
include("pt_scalars_pre.jl")

include("induction.jl")
include("lorentz.jl")
include("lorentz_pre.jl")
include("induction_pre.jl")


end