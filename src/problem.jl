export LimaceProblem

import Base: iterate, length

abstract type Forcing{T} end

length(f::Forcing) = 1
iterate(f::Forcing) = (f, nothing)
iterate(f::Forcing, ::Any) = nothing

Base.show(io::IO, f::T) where T <: Forcing = print(io, "$T(N = $(getfield(f,1).N), factor = $(f.factor))")

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct Coriolis{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct Inertial{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct Diffusion{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct InductionU0{TB,T} <: Forcing{1}
    basis::Basis{TB}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct InductionB0{TB,TU,T} <: Forcing{2}
    bbasis::Basis{TB}
    ubasis::Basis{TU}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct Lorentz{TU,TB,T} <: Forcing{2}
    ubasis::Basis{TU}
    bbasis::Basis{TB}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

"""
$(TYPEDEF)

$(TYPEDFIELDS)
"""
mutable struct Advection{TU,T} <: Forcing{1}
    basis::Basis{TU}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function Inertial(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Inertial(b, ComplexF64(factor), mat, false)
end

function Coriolis(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Coriolis(b, ComplexF64(factor), mat, false)
end

function Diffusion(b::Basis, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return Diffusion(b, ComplexF64(factor), mat, false)
end

function InductionU0(b::Basis, U0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(b), length(b))
    return InductionU0(b, U0, ComplexF64(factor), mat, false)
end

function InductionB0(bb::Basis, ub::Basis, B0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(bb), length(ub))
    return InductionB0(bb, ub, B0, ComplexF64(factor), mat, false)
end

function Lorentz(ub::Basis, bb::Basis, B0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(ub), length(bb))
    return Lorentz(ub, bb, B0, ComplexF64(factor), mat, false)
end

function Advection(ub::Basis, U0, factor::T=1.0) where T
    mat = spzeros(ComplexF64, length(ub), length(ub))
    return Advection(ub, U0, ComplexF64(factor), mat, false)
end


function assemble!(f::Inertial; kwargs...)
    f.mat = sparse(Limace.inertial(f.basis; kwargs...))
    f.preassembled = true
    return f.mat
end

function assemble!(f::Coriolis; kwargs...)
    f.mat = Limace.coriolis(f.basis; kwargs...)
    f.preassembled = true
    return f.mat
end

function assemble!(f::Diffusion; kwargs...)
    f.mat = Limace.diffusion(f.basis; kwargs...)
    f.preassembled = true
    return f.mat
end

function assemble!(f::InductionU0; kwargs...)
    f.mat = sum(Limace.induction(f.basis, U0, f.basis; kwargs...) for U0 in f.U0)
    f.preassembled = true
    return f.mat
end

function assemble!(f::InductionB0; kwargs...)
    f.mat = sum(Limace.induction(f.bbasis, f.ubasis, B0; kwargs...) for B0 in f.B0)
    f.preassembled = true
    return f.mat
end

function assemble!(f::Lorentz; kwargs...)
    f.mat = sum(Limace.lorentz(f.ubasis, f.bbasis, B0; kwargs...) for B0 in f.B0)
    f.preassembled = true
    return f.mat
end

function assemble!(f::Advection; kwargs...)
    f.mat = sum(Limace.advection(f.basis, U0; kwargs...) for U0 in f.U0)
    f.preassembled = true
    return f.mat
end


## LimaceProblem

"""
$(TYPEDEF)

$(TYPEDFIELDS)

"""
mutable struct LimaceProblem{T}
    bases
    forcings 
    RHS::SparseMatrixCSC{T}
    LHS::SparseMatrixCSC{T}
    sol::Union{Eigen{T, T, Matrix{T}, Vector{T}}, GeneralizedEigen{T, T, Matrix{T}, Vector{T}}}
    preassembled::Bool
    assembled::Bool
    solved::Bool
end


function LimaceProblem(bases, forcings=Forcing[])
    RHS = spzeros(ComplexF64, sum(length.(bases)), sum(length.(bases)))
    LHS = spzeros(ComplexF64, sum(length.(bases)), sum(length.(bases)))
    sol = Eigen(ComplexF64[], Matrix{ComplexF64}(undef, 0 , 0))
    return LimaceProblem{ComplexF64}(bases, forcings, RHS, LHS, sol, false, false, false)
end

# function Base.show(io::IO, problem::LimaceProblem) 
#     print(io, "LimaceProblem(bases = $([typeof(b) for b in problem.bases]), forcings = $(problem.forcings), preassembled = $(problem.preassembled), assembled = $(problem.assembled), solved = $(problem.solved))")
# end
## assemble the problem

"""
    preassemble!(problem::LimaceProblem; kwargs...)

Preassemble `problem.forcing` matrices in the problem.
"""
function preassemble!(problem::LimaceProblem; kwargs...)
    for f in problem.forcings
        if !f.preassembled
            assemble!(f; kwargs...)
        end
    end
    problem.preassembled=true
    return nothing
end

"""
    assemble!(problem::LimaceProblem; kwargs...)

Assemble the problem matrices `problem.LHS` and `problem.RHS` from the forcing matrices that may or may not be preassembled.
For now, only `Limace.Inertial` are added to the `LHS` matrix.
"""
function assemble!(problem::LimaceProblem; kwargs...)
    if !problem.preassembled 
        preassemble!(problem; kwargs...)
    end

    prematLHS = [spzeros(ComplexF64,length(bi),length(bj)) for bi in problem.bases, bj in problem.bases]
    prematRHS = [spzeros(ComplexF64,length(bi),length(bj)) for bi in problem.bases, bj in problem.bases]

    for f in problem.forcings
        if typeof(f)<:Limace.Inertial
            _add_to_premat!(problem, prematLHS, f)
        else
            _add_to_premat!(problem, prematRHS, f)
        end
    end

    problem.LHS = hvcat(length(problem.bases),permutedims(prematLHS)...)
    problem.RHS = hvcat(length(problem.bases),permutedims(prematRHS)...)
    problem.assembled=true
    return problem.LHS, problem.RHS
end

"""
    add!(problem::LimaceProblem, f::Forcing)

Add a forcing `f` to the `problem`. The forcing is not preassembled.
"""
function add!(problem::LimaceProblem, f::Forcing)
    push!(problem.forcings, f)
    problem.preassembled = false
    return nothing
end

function _add_to_premat!(problem::LimaceProblem, premat, f::TF) where {TF <: Forcing{1}}
    basis = getfield(f,1)
    ib = last(findfirst(isequal(basis),problem.bases))
    premat[ib,ib] += f.mat*f.factor
    return nothing
end

function _add_to_premat!(problem::LimaceProblem, premat, f::TF) where {TF <: Forcing{2}}
    basis1 = getfield(f,1)
    basis2 = getfield(f,2)
    ib1 = last(findfirst(isequal(basis1),problem.bases))
    ib2 = last(findfirst(isequal(basis2),problem.bases))
    premat[ib1,ib2] += f.mat*f.factor
    return nothing
end


## Solving the problem

"""
    solve!(problem::LimaceProblem; method=:dense, kwargs...)

Solve `problem` using the specified method. The default method is `:dense`, which transforms the problem matrices to dense matrices.
Other methods are `:sparse`, which uses the sparse matrices directly.

Solutions are stored in `problem.sol` and `problem.solved` is set to `true`.
"""
function solve!(problem::LimaceProblem; method=:dense, kwargs...)
    if !problem.assembled
        @warn "Problem was not assembled. Assembling the problem now."
        assemble!(problem)
    end

    if method == :dense
        return solve_dense!(problem)
    elseif method == :sparse
        return solve_sparse!(problem; kwargs...)
    else
        error("Unknown method $(method). Use :dense or :sparse.")
    end
  
    return problem.sol
end


function solve_dense!(problem::LimaceProblem)

    if (first(problem.bases).N > 30) && length(first(problem.bases).m) > 1
        @warn "LimaceProblem.solve_dense! is not optimized for large problems. Use :sparse method instead."
    end

    if isdiag(problem.LHS)
        if problem.LHS ≈ I
            C = Matrix(problem.RHS)
        else
            C = Matrix(problem.RHS\Diagonal(problem.LHS))
        end
        problem.sol = eigen(C)
    else
        problem.sol = eigen(Matrix(problem.RHS), Matrix(problem.LHS))
    end
    problem.solved = true

    return problem.sol
end

function solve_sparse!(problem::LimaceProblem; target=Inf, kwargs...)
    if isinf(target) 
        λ, x = EigenSolve.eigs(problem.RHS, problem.LHS; kwargs...)
    else
        λ, x = EigenSolve.eigstarget(problem.RHS, problem.LHS, target; kwargs...)
    end

    problem.sol = GeneralizedEigen(λ, x)
    problem.solved = true

    return problem.sol
end