module Limace

using SparseArrays
using LinearAlgebra
using DocStringExtensions
using Reexport

include("Utils.jl")
using .Utils

include("Poly.jl")
using .Poly

include("Quad.jl")
using .Quadrature

include("Bases.jl")
@reexport using .Bases
using .Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax

# forces
include("forces/inertial.jl")
include("forces/coriolis.jl")
include("forces/diffusion.jl")
include("forces/induction.jl")
include("forces/lorentz.jl")

include("bases/Inviscid.jl")
@reexport using .InviscidBasis

include("bases/inviscid2.jl")
using .InviscidBasis2

include("bases/InviscidNoBC.jl")
using .InviscidBasisNoBC

include("bases/Insulating.jl")
@reexport using .InsulatingBasis

include("forces/specializations.jl")
# include("bases/viscous.jl")
# include("bases/viscous_chen.jl")

# include("discrete/DiscretePart.jl")

# include("mhdproblem.jl")

end
