module EigenProblem

using ArnoldiMethod, LinearMaps, LinearAlgebra, SparseArrays

export eigstarget

#Use Shift-Invert method to solve generalized eigen problem, factorization is done using UMFPACK (lu).
#λBx = Ax -> 1/(λ-σ)x = (A-σB)⁻¹Bx
function eigstarget(A, B, σ; kwargs...)
    C = A - σ * B
    P = lu(C)
    LO = LinearMap{eltype(C)}((y, x) -> ldiv!(y, P, B * x), size(C, 2))
    pschur, history = partialschur(LO; kwargs...)
    evals, x = partialeigen(pschur)
    λ = 1 ./ evals .+ σ
    return λ, x
end



end #module