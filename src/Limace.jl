module Limace

using SparseArrays
using LinearAlgebra
using DocStringExtensions
using Reexport

include("Utils.jl")
using .Utils

include("Poly.jl")
using .Poly

include("Bases.jl")
@reexport using .Bases
using .Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p_l, lmn_t_l, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax


include("Quad.jl")
using .Quadrature

include("Discretization.jl")
using .Discretization

# forces
include("forces/inertial.jl")
include("forces/coriolis.jl")
include("forces/diffusion.jl")
include("forces/induction.jl")
include("forces/lorentz.jl")

include("bases/Inviscid.jl")
@reexport using .InviscidBasis

include("bases/InviscidShell.jl")
using .InviscidShellBasis

include("bases/Unconstrained.jl")
using .UnconstrainedBasis

include("bases/InviscidNoBC.jl")
using .InviscidBasisNoBC

include("bases/Insulating.jl")
@reexport using .InsulatingBasis

include("bases/InsulatingNoBC.jl")
using .InsulatingBasisNoBC

include("bases/Viscous.jl")
@reexport using .ViscousBasis

include("bases/ViscousNoBC.jl")
@reexport using .ViscousBasisNoBC

include("bases/ThinWall.jl")
using .ThinWallBC

include("bases/ViscousShell.jl")
@reexport using .ViscousShellBasis

include("forces/specializations.jl")


end
