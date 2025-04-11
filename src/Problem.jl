export LimaceProblem



import Base: iterate, length

abstract type Forcing{T} end

length(f::Forcing) = 1
iterate(f::Forcing) = (f, nothing)
iterate(f::Forcing, ::Any) = nothing

Base.show(io::IO, f::T) where T <: Forcing = print(io, "$T(N = $(getfield(f,1).N), factor = $(f.factor))")

mutable struct Coriolis{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

mutable struct Inertial{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

mutable struct Diffusion{TB,T} <: Forcing{1}
    basis::Basis{TB}
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

mutable struct InductionU0{TB,T} <: Forcing{1}
    basis::Basis{TB}
    U0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

mutable struct InductionB0{TB,TU,T} <: Forcing{2}
    bbasis::Basis{TB}
    ubasis::Basis{TU}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

mutable struct Lorentz{TU,TB,T} <: Forcing{2}
    ubasis::Basis{TU}
    bbasis::Basis{TB}
    B0
    factor::T
    mat::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

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


mutable struct LimaceProblem 
    bases::AbstractVector{Basis}
    forcings # [ [[f1_b1_b1, f2_b1_b1], f_b1_b2, ...] , [f_b2_b2, [f1_b2_b1, f2_b2_b3]] , ... ]
    RHS::SparseMatrixCSC{ComplexF64}
    LHS::SparseMatrixCSC{ComplexF64}
    preassembled::Bool
end

function LimaceProblem(bases, forcings=Forcing[])
    RHS = spzeros(ComplexF64, sum(length.(bases)), sum(length.(bases)))
    LHS = spzeros(ComplexF64, sum(length.(bases)), sum(length.(bases)))
    return LimaceProblem(bases, forcings, RHS, LHS, false)
end

function preassemble!(problem::LimaceProblem)
    for f in problem.forcings
        if !f.preassembled
            assemble!(f)
        end
    end
    problem.preassembled=true
    return nothing
end

function assemble!(problem::LimaceProblem)
    if !problem.preassembled 
        preassemble!(problem)
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
    return problem.LHS, problem.RHS
end

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




