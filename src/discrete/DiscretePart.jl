module DiscretePart

using ForwardDiff
using SparseArrays
using Wigxjpf
using FastGaussQuadrature

#only for dev
wig_table_init(200, 9)
wig_temp_init(200)

include("poly.jl")
include("quad.jl")
include("pt_scalars.jl")

include("induction.jl")


end