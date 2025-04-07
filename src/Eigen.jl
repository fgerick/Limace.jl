module Eigen

using DocStringExtensions

using ArnoldiMethod, LinearMaps, LinearAlgebra, SparseArrays


"""
$(TYPEDSIGNATURES)

Use Shift-Invert method to solve generalized eigen problem, factorization is done using UMFPACK (lu).

```math
\\lambda\\mathbf{B}\\mathbf{x} = \\mathbf{A}\\mathbf{x} \\rightarrow \\frac{1}{\\lambda-\\sigma}\\mathbf{x} = (\\mathbf{A}-\\sigma\\mathbf{B})^{-1}\\mathbf{B}\\mathbf{x}.
```

Keyword arguments are passed to `partialschur` function (see [ArnoldiMethod.partialschur](https://julialinearalgebra.github.io/ArnoldiMethod.jl/stable/#ArnoldiMethod.partialschur) for details).
The function returns the eigenvalues `λ` and eigenvectors `x` of the generalized eigenvalue problem.
"""
function eigstarget(A, B, σ; kwargs...)
    C = A - σ * B
    P = lu(C)
    LO = LinearMap{eltype(C)}((y, x) -> ldiv!(y, P, B * x), size(C, 2))
    pschur, history = partialschur(LO; kwargs...)
    evals, x = partialeigen(pschur)
    λ = 1 ./ evals .+ σ
    return λ, x
end

"""
$(TYPEDSIGNATURES)

Use Inverse method to solve generalized eigen problem, factorization is done using UMFPACK (lu).

```math
\\lambda\\mathbf{B}\\mathbf{x} = \\mathbf{A}\\mathbf{x} \\rightarrow \\lambda\\mathbf{x} = \\mathbf{B}^{-1}\\mathbf{A}\\mathbf{x}.
```

Keyword arguments are passed to `partialschur` function (see [ArnoldiMethod.partialschur](https://julialinearalgebra.github.io/ArnoldiMethod.jl/stable/#ArnoldiMethod.partialschur) for details).
The function returns the eigenvalues `λ` and eigenvectors `x` of the generalized eigenvalue problem.
"""
function eigs(A, B; kwargs...)
    P = lu(B)
    LO = LinearMap{eltype(A)}((y, x) -> ldiv!(y, P, A * x), size(A, 2))
    pschur, history = partialschur(LO; kwargs...)
    λ, x = partialeigen(pschur)
    return λ, x
end

function eigs(A; kwargs...)
    pschur, history = partialschur(A; kwargs...)
    λ, x = partialeigen(pschur)
    return λ, x
end


end #module