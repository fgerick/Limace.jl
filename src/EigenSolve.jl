module EigenSolve

using DocStringExtensions

using ArnoldiMethod
using LinearMaps

using LinearAlgebra
using SparseArrays
using Statistics

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


# (pslast, RHS, LHS, target; nev)
function _eigstargetumfpack(P,B, σ; kwargs...)
    LO = LinearMap{ComplexF64}((y,x)->ldiv!(y,P,B*x),size(B,2))
    pschur, _ = partialschur(LO; kwargs...)
    evals, x = partialeigen(pschur)
    λ = 1 ./evals .+ σ 
    return λ,x
end

function _eigstargetumfpack(A,B,C,ps,σ; kwargs...)
	for (i,j,_) in zip(findnz(C)...)
		C[i,j] = A[i,j]-σ*B[i,j]
	end
	lu!(ps, C; check=false, reuse_symbolic=true)
    LO = LinearMap{ComplexF64}((y,x)->ldiv!(y,ps,B*x),size(A,2))
    pschur, _ = partialschur(LO; kwargs...)
    evals, x = partialeigen(pschur)
    λ = 1 ./evals .+ σ 
    return λ,x
end

"""
$(TYPEDSIGNATURES)

Track eigenvalues `targets0` and corresponding eigenvectors `evecs0`
"""
function tracking(targets0::Vector{ComplexF64}, evecs0, LHS, RHS0, updateRHS!, α0, αmax; 
	stepthresh = sqrt(eps()), 
	last_special=false, 
	nev=5, 
	corthresh=0.95, 
	δα0 = (αmax-α0)/20, 
	maxstep=abs(αmax-α0)/10,
	itermax=10_000, 
	info=false
	)
	    
	direction_positive = δα0 > 0
	if direction_positive
		check_arrival = >=
	else
		check_arrival = <=
	end

	evals = [ComplexF64[] for target in targets0]
	evecs = [Vector{ComplexF64}[] for u in evecs0]
	corrs = [Float64[] for _ in 1:length(targets0)]
	
	αs = [Float64[] for _ in 1:length(targets0)]
	
	
	RHS = copy(RHS0)
	_C = RHS - randn(ComplexF64)*LHS #random eigenvalue
	ps = lu(_C)

	for (itarget,(target0, u0)) in enumerate(zip(targets0,evecs0))
		
		@info "target $(itarget)/$(length(targets0))"

		α = α0
		αt = α0
		δα = δα0
		iter = 0

		push!(evals[itarget],target0)
		push!(evecs[itarget],u0)
		push!(αs[itarget], α0)
		target = target0

		while !check_arrival(α,αmax) && (iter <= itermax)
			αt = α + δα
			if check_arrival(αt,αmax)
				αt = αmax
			end

			if last_special && (αt ≈ αmax)
				updateRHS!(RHS, RHS0, αt, α0)
				pslast = lu(RHS-target*LHS)
				evals1, evecs1 = _eigstargetumfpack(pslast,LHS, target; nev)
			else
				updateRHS!(RHS, RHS0, αt, α0)	
				evals1, evecs1 = try
					lu!(ps, RHS-target*LHS)
					_eigstargetumfpack(RHS,LHS,_C, ps, target; nev)
				catch e
					@warn "error eigensolution, next target"
					throw(e)
					break
				end
			end

			
			_corrs = [abs(cor(@views(evecs1[:,i]),evecs[itarget][end])) for i in axes(evecs1,2)]
			corrmax, imax = findmax(_corrs)
			info && @info "α = $(αt), max. correlation = $corrmax"
			if (corrmax <= corthresh)
				δα/=2
				info && @info "δα = $(δα)"
				if abs(δα) <= stepthresh
					@warn "stuck, stopping"
					break
				end
				continue
			else
				push!(corrs[itarget],corrmax)
				push!(evals[itarget], evals1[imax])
				push!(evecs[itarget], evecs1[:,imax])
				push!(αs[itarget],αt)
				target = evals1[imax]
				α = αt
				if abs(δα) < maxstep
					δα*=1.5
				end
			end
			iter +=1
		end


	end
	return corrs, evals, evecs, αs
end 

end #module