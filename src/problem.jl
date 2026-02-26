export LimaceProblem

import Base: iterate, length

abstract type Forcing{T} end

length(f::Forcing) = 1
iterate(f::Forcing) = (f, nothing)
iterate(f::Forcing, ::Any) = nothing

Base.show(io::IO, f::T) where T <: Forcing = print(io, "$T(N = $(getfield(f,1).N), factor = $(f.factor))")

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
    sol::Union{Eigen{T}, GeneralizedEigen{T}}
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
    preassemble!(problem::LimaceProblem; threads=false, kwargs...)

Preassemble `problem.forcing` matrices in the problem. When `threads=true` the assembly is done using `Threads.nthreads()` threads.
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
    assemble!(problem::LimaceProblem; threads=false, kwargs...)

Assemble the problem matrices `problem.LHS` and `problem.RHS` from the forcing matrices that may or may not be preassembled.
For now, only `Limace.Inertial` are added to the `LHS` matrix. When `threads=true` the assembly is done using `Threads.nthreads()` threads.
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
    if abs(f.factor) != 0
        basis = getfield(f,1)
        ib = last(findfirst(isequal(basis),problem.bases))
        premat[ib,ib] += f.mat*f.factor
    end
    return nothing
end

function _add_to_premat!(problem::LimaceProblem, premat, f::TF) where {TF <: Forcing{2}}
    if abs(f.factor) != 0
        basis1 = getfield(f,1)
        basis2 = getfield(f,2)
        ib1 = last(findfirst(isequal(basis1),problem.bases))
        ib2 = last(findfirst(isequal(basis2),problem.bases))
        premat[ib1,ib2] += f.mat*f.factor
    end
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
            C = Matrix(Diagonal(problem.LHS)\problem.RHS)
        end
        problem.sol = eigen(C)
    else
        problem.sol = eigen(Matrix(problem.RHS), Matrix(problem.LHS))
    end
    problem.solved = true

    return problem.sol
end

function solve_sparse!(problem::LimaceProblem; target=Inf, kwargs...)
    if typeof(target) <: Number
        if isinf(target) 
            λ, x = EigenSolve.eigs(problem.RHS, problem.LHS; kwargs...)
        else
            λ, x = EigenSolve.eigstarget(problem.RHS, problem.LHS, target; kwargs...)
        end
    elseif typeof(target)<:AbstractVector
        Tc = complex(eltype(problem.LHS))
        C = problem.RHS - first(target)*problem.LHS
        P = lu(C)
        λ, x = EigenSolve._eigstargetumfpack(P, problem.LHS, first(target); kwargs...)
        for t in target[2:end]
            _λ, _x = EigenSolve._eigstargetumfpack(problem.RHS, problem.LHS, C, P, t; kwargs...)
            for (i,λi) in enumerate(_λ)
                if !any(isapprox(λi, atol=10sqrt(eps())), λ)
                    append!(λ,λi)
                    @views x = hcat(x,_x[:,i])
                end
            end
        end
    else
        @error "target should be a complex or real number, or an AbstractVector of real/complex numbers."
    end

    problem.sol = GeneralizedEigen(λ, x)
    problem.solved = true

    return problem.sol
end