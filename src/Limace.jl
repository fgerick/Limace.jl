module Limace

using LinearAlgebra
using DocStringExtensions
using SparseArrays
using Reexport

include("Utils.jl")
using .Utils

include("Poly.jl")
using .Poly

include("Bases.jl")
@reexport using .Bases
using .Bases: nrange_p, nrange_t, nrange_p_bc, nrange_t_bc, np, nt, t, s, bcs_p, bcs_t, lmn_p, lmn_t, lmn2k_p_dict, lmn2k_t_dict, lpmax, ltmax

include("Quadrature.jl")
using .Quadrature

include("Discretization.jl")
using .Discretization

include("EigenSolve.jl")
using .EigenSolve

include("problem.jl")

# forces
include("forces/advection.jl")
include("forces/inertial.jl")
include("forces/coriolis.jl")
include("forces/diffusion.jl")
include("forces/crossterm.jl")
include("forces/induction.jl")
include("forces/lorentz.jl")
include("forces/scalaradvection.jl")
include("forces/buoyancy.jl")
include("forces/bc.jl")

# bases

include("bases/Inviscid.jl")
@reexport using .InviscidBasis

include("bases/PerfectlyConducting.jl")
@reexport using .PerfectlyConductingBasis

include("bases/Insulating.jl")
@reexport using .InsulatingBasis

include("bases/Viscous.jl")
@reexport using .ViscousBasis

include("bases/InviscidShell.jl")
using .InviscidShellBasis

include("bases/Unconstrained.jl")
@reexport using .UnconstrainedBasis

include("bases/InviscidNoBC.jl")
using .InviscidBasisNoBC

include("bases/InsulatingNoBC.jl")
using .InsulatingBasisNoBC

include("bases/ViscousNoBC.jl")
using .ViscousBasisNoBC

include("bases/ThinWall.jl")
using .ThinWallBasis

include("bases/ViscousShell.jl")
using .ViscousShellBasis

include("bases/Temperature.jl")
using .TemperatureBasis

include("forces/specializations.jl")

include("Processing.jl")
using .Processing

end
