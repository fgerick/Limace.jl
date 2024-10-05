module Processing

using DocStringExtensions

using LinearAlgebra
using SparseArrays
using Statistics
using ..Bases
using ..Bases: lmn_p, lmn_t, lmn_p_bc, lmn_t_bc, lmn2k_p_dict, lmn2k_t_dict, _lmn2cdeg_p, _lmn2cdeg_t



function getindices(u0, b0, u1, b1)
    lmnp0 = lmn_p(u0)
    lmnt0 = lmn_t(u0)
    np0 = length(lmnp0)
    nu = np0 + length(lmnt0)
    lmnpnew = lmn_p(u1)
    lmntnew = lmn_t(u1)
    indsp = findall(in(lmnpnew), lmnp0)
    indst = findall(in(lmntnew), lmnt0) .+ np0
    indsu = vcat(indsp, indst)

    lmnpnew = lmn_p_bc(u1)
    lmntnew = lmn_t_bc(u1)
    indsp = findall(in(lmnpnew), lmnp0)
    indst = findall(in(lmntnew), lmnt0) .+ np0
    indsu_bc = vcat(indsp, indst)

    lmnp0 = lmn_p(b0)
    lmnt0 = lmn_t(b0)
    np0 = length(lmnp0)

    lmnpnew = lmn_p(b1)
    lmntnew = lmn_t(b1)
    indsp = findall(in(lmnpnew), lmnp0) .+ nu
    indst = findall(in(lmntnew), lmnt0) .+ np0 .+ nu
    indsb = vcat(indsp, indst)

    lmnpnew = lmn_p_bc(b1)
    lmntnew = lmn_t_bc(b1)
    indsp = findall(in(lmnpnew), lmnp0) .+ nu
    indst = findall(in(lmntnew), lmnt0) .+ np0 .+ nu
    indsb_bc = vcat(indsp, indst)



    return vcat(indsu_bc, indsb_bc), vcat(indsu, indsb)
end


function getindices(u0, u1)
    lmnp0 = lmn_p(u0)
    lmnt0 = lmn_t(u0)
    np0 = length(lmnp0)
    lmnpnew = lmn_p(u1)
    lmntnew = lmn_t(u1)
    indsp = findall(in(lmnpnew), lmnp0)
    indst = findall(in(lmntnew), lmnt0) .+ np0
    indsu = vcat(indsp, indst)

    lmnpnew = lmn_p_bc(u1)
    lmntnew = lmn_t_bc(u1)
    indsp = findall(in(lmnpnew), lmnp0)
    indst = findall(in(lmntnew), lmnt0) .+ np0
    indsu_bc = vcat(indsp, indst)

    return indsu_bc, indsu
end

function getindices_nobc(u, b)

    lmnpu = lmn_p_bc(u)
    lmnpud = lmn2k_p_dict(u)
    ipu = [lmnpud[i] for i in lmnpu]
    jpu = [lmnpud[i] for i in lmn_p(u)]
    np0 = length(lmn_p(u))

    lmntu = lmn_t_bc(u)
    lmntud = lmn2k_t_dict(u)
    itu = [lmntud[i] + np0 for i in lmntu]
    iu = vcat(ipu, itu)

    jtu = [lmntud[i] + np0 for i in lmn_t(u)]
    ju = vcat(jpu, jtu)

    nu = np0 + length(lmn_t(u))

    lmnpb = lmn_p_bc(b)
    lmnpbd = lmn2k_p_dict(b)
    ipb = [lmnpbd[i] + nu for i in lmnpb]
    jpb = [lmnpbd[i] + nu for i in lmn_p(b)]
    lmntb = lmn_t_bc(b)
    lmntbd = lmn2k_t_dict(b)

    np0 = length(lmn_p(b))
    itb = [lmntbd[i] + np0 + nu for i in lmntb]
    jtb = [lmntbd[i] + np0 + nu for i in lmn_t(b)]
    ib = vcat(ipb, itb)
    jb = vcat(jpb, jtb)

    return vcat(iu, ib), vcat(ju, jb)
end

function reducematrix(a, u0, b0, u1, b1)
    i, j = getindices(u0, b0, u1, b1)
    n = length(u1) + length(b1)
    out = spzeros(ComplexF64, n, n)
    a1 = a[i, j]
    i2, j2 = getindices_nobc(u1, b1)
    out[i2, j2] = a1
    return out
end

function usection(::Val{true}, λ, λ2, evecs, evecs2, u0, b0, u1, b1; threshc=0.99, threshλ=0.1)
    @assert u1.N > u0.N
    nth = Threads.nthreads()
    is = [Int[] for _ in 1:nth]
    is2 = [Int[] for _ in 1:nth]
    inds = getindices(u1, b1, u0, b0)[2]

    Threads.@threads for i in eachindex(λ)
        it = Threads.threadid()
        found = false
        ω = λ[i]
        for i2 in eachindex(λ2)
            ω2 = λ2[i2]
            if (!found) && (abs(ω) * (1 - threshλ) < abs(ω2) < (1 + threshλ) * abs(ω))
                c = cor(@views(evecs[:, i]), @views(evecs2[inds, i2]))
                if (abs(c) > threshc) && (sign(imag(ω)) == sign(imag(ω2)))
                    push!(is[it], i)
                    push!(is2[it], i2)
                    found = true
                    @goto λi
                end
            end
        end
        @label λi
    end
    return vcat(is...), vcat(is2...)
end

"""
$(TYPEDSIGNATURES)

Compute the kinetic and magnetic energies of all eigenvectors `us`, where `size(us,2)` is the number of eigenvalues. 
`LHS` is the mass-matrix and `u` the velocity basis used in the calculation of the eigenvectors `us`.
"""
function ekinmags(us, LHS, u)
    nu = length(u)
    nm = size(LHS, 1)
    LHSu = Diagonal(view(LHS, 1:nu, 1:nu))
    LHSb = SymTridiagonal(view(LHS, nu+1:nm, nu+1:nm))
    nev = size(us, 2)
    ekin = zeros(nev)
    emag = zeros(nev)
    Threads.@threads for j in axes(us, 2)
        ekin[j] = abs(dot(view(us, 1:nu, j), LHSu, view(us, 1:nu, j)))
        emag[j] = abs(dot(view(us, nu+1:nm, j), LHSb, view(us, nu+1:nm, j)))
    end
    return ekin, emag
end

"""
$(TYPEDSIGNATURES)

Calculate spectrum for all eigenvectors `evecs` comprised of velocity basis `u` and magnetic field basis `b`. 
Use keyword `lmn=1,2,3` to select `l=1`, `m=2` or `n=3`.
"""
function spectrum(evecs, u, b; lmn=1)
    N = max(u.N, b.N)
    _N = N + 1
    nev = size(evecs, 2)
    spec = zeros(4 * _N, nev)
    spec_fac = zeros(Int, 4 * _N)


    lmnpu = lmn_p(u)
    np = length(lmnpu)
    lmntu = lmn_t(u)
    nt = length(lmntu)
    lmnpb = lmn_p(b)
    npb = length(lmnpb)
    lmntb = lmn_t(b)

    nu = np + nt
    j = 1
    for k in axes(evecs, 1)
        if k <= np
            j = abs(lmnpu[k][lmn]) + 1
        elseif k <= nu
            j = abs(lmntu[k-np][lmn]) + 1 + _N
        elseif k <= nu + npb
            j = abs(lmnpb[k-nu][lmn]) + 1 + 2 * _N
        else
            j = abs(lmntb[k-nu-npb][lmn]) + 1 + 3 * _N
        end
        for i in axes(evecs, 2)
            α = abs(evecs[k, i])^2
            spec[j, i] += α
        end
        spec_fac[j] += 1
    end

    spec_fac[spec_fac.==0] .= 1
    spec ./= spec_fac

    specup, specut, specbp, specbt = spec[1:_N, :], spec[_N+1:2*_N, :], spec[2*_N+1:3*_N, :], spec[3*_N+1:4*_N, :]
    return specup, specut, specbp, specbt
end

# function spec_ub_all_lmn(evecs, u, b; lmn=1)
#     N = max(u.N, b.N)
#     _N = N + 1
#     nev = size(evecs, 2)

#     specup, specut = zeros(nev,_N), zeros(nev,_N)
#     specbp, specbt = zeros(nev,_N), zeros(nev,_N)

#     lmnpu = lmn_p(u)
#     np = length(lmnpu)
#     lmntu = lmn_t(u)
#     nt = length(lmntu)
#     lmnpb = lmn_p(b)
#     npb = length(lmnpb)
#     lmntb = lmn_t(b)	
#     nu = np + nt

#     for i in axes(evecs,2)
#         k = 1
#         for k in axes(evecs,1)
#             α = abs(evecs[k,i])^2
#             if k<=np
#                 specup[i,abs(lmnpu[k][lmn])+1]+=α
#             elseif k<=nu
#                 specut[i,abs(lmntu[k-np][lmn])+1]+=α
#             elseif k<=nu+npb
#                 specbp[i,abs(lmnpb[k-nu][lmn])+1]+=α
#             else
#                 specbt[i,abs(lmntb[k-nu-npb][lmn])+1]+=α
#             end
#         end
#     end

#     return specup, specut, specbp, specbt
# end

"""
$(TYPEDSIGNATURES)

	Get all `(l,m,n)` that correspond to the poloidal and toroidal component of `u` and `b` basis at each Cartesian
	degree `ñ ∈ 1:N`.
"""
function lmn_n(u::Basis{TU}, b::Basis{TB}) where {TU, TB}
    N = u.N
    Ns = 1:N
    _u(N) = TU(N; m=u.m)
    _b(N) = TB(N; m=b.m)
    lmnpu = [setdiff(lmn_p(_u(N)), lmn_p(_u(N - 1))) for N in Ns]
    lmntu = [setdiff(lmn_t(_u(N)), lmn_t(_u(N - 1))) for N in Ns]
    lmnpb = [setdiff(lmn_p(_b(N)), lmn_p(_b(N - 1))) for N in Ns]
    lmntb = [setdiff(lmn_t(_b(N)), lmn_t(_b(N - 1))) for N in Ns]

    return Ns, lmnpu, lmntu, lmnpb, lmntb
end


"""
$(TYPEDSIGNATURES)

Compute poloidal/toroidal kinetic/magnetic spectra of `evecs` as a function of max. Cartesian monomial degree.
`evecs` are the eigenvectors computed using the velocity basis `u` and magnetic field basis `b`.
"""
function spectrum_cartesian(evecs, u, b)
    N = max(u.N, b.N)
    nev = size(evecs, 2)
    spec = zeros(4N, nev)
    spec_fac = zeros(Int, 4N)

    lmnpu = lmn_p(u)
    np = length(lmnpu)
    lmntu = lmn_t(u)
    nt = length(lmntu)
    lmnpb = lmn_p(b)
    npb = length(lmnpb)
    lmntb = lmn_t(b)

    nu = np + nt
    j = 1
    @inbounds for k in axes(evecs, 1)
        if k <= np
            j = _lmn2cdeg_p(u,lmnpu[k]...) 
        elseif k <= nu
            j = _lmn2cdeg_t(u,lmntu[k-np]...) + N 
        elseif k <= nu + npb
            j = _lmn2cdeg_p(b,lmnpb[k-nu]...) + 2N 
        else
            j = _lmn2cdeg_t(b,lmntb[k-nu-npb]...) + 3N 
        end
        for i in axes(evecs, 2)
            α = abs(evecs[k, i])^2
            spec[j, i] += α
        end
        spec_fac[j] += 1
    end
    spec_fac[spec_fac.==0] .= 1
    spec ./= spec_fac

    specup, specut, specbp, specbt = spec[1:N, :], spec[N+1:2N, :], spec[2N+1:3N, :], spec[3N+1:4N, :]
    return specup, specut, specbp, specbt
end

"""
$(TYPEDSIGNATURES)

Compute the ratio of peak energy to energy at truncation degree
(max between two last Cartesian degrees) in toroidal/poloidal kinetic/magnetic energy.
`evecs` are the eigenvectors computed using the velocity basis `u` and magnetic field basis `b`.
"""
function epeak_etrunc_cartesian(evecs, u, b)
    ratios = zeros(4, size(evecs, 2))

    specs = spectrum_cartesian(evecs, u, b)
    maxima = maximum(vcat(maximum.(specs, dims=1)...), dims=1)[:]
    for (i, spec) in enumerate(specs)
        @views ratios[i, :] .= (maximum(spec[end-1:end, :], dims=1)[:] ./ maxima)
    end
    return ratios
end

"""
$(TYPEDSIGNATURES)

Find all indices of `evals1` and `evals2` for which `findall(y->any(x->isapprox(x,y; rtol=λtol), evals2),evals1)`.
Multithreaded.
"""
function eigenvalue_filter(evals1, evals2; λtol=1e-3)
    nt = Threads.nthreads()
    is1 = [Int[] for _ in 1:nt]
    is2 = [Int[] for _ in 1:nt]

    Threads.@threads for i1 in eachindex(evals1)
        λ1 = evals1[i1]
        δ, i2 = findmin(x -> abs(x - λ1), evals2)
        if δ / abs(λ1) < λtol
            id = Threads.threadid()
            push!(is1[id], i1)
            push!(is2[id], i2)
        end
    end
    return vcat(is1...), vcat(is2...)
end

"""
$(TYPEDSIGNATURES)

Find numerically converged eigensolutions between two resolutions.
"""
function numerical_filter(evecs1, evecs2, evals1, evals2, u1, u2, b1, b2;
    threshc=0.95,
    λtol=1e-2,
    threads=true,
    epeakratiothresh=1e-1,
    kwargs...
	)

	@assert size(evecs1,1) == length(u1)+length(b1)
	@assert size(evecs2,1) == length(u2)+length(b2)

    is1, is2 = eigenvalue_filter(evals1, evals2; λtol)

    ratios1 = maximum(epeak_etrunc_cartesian(@views(evecs1[:, is1]), u1, b1), dims=1)[:]
    _s11 = ratios1 .< epeakratiothresh
    ratios2 = maximum(epeak_etrunc_cartesian(@views(evecs2[:, is2]), u2, b2), dims=1)[:]
    _s22 = ratios2 .< epeakratiothresh

    _is1 = is1[_s11]
    _is2 = is2[_s22]

	if threshc > 0.0
		tval = Val(threads)
		@views is1, is2 = usection(tval, evals1[_is1], evals2[_is2], evecs1[:, _is1], evecs2[:, _is2], u1, b1, u2, b2; threshc, kwargs...)

		_is1 = _is1[is1]
		_is2 = _is2[is2]
	end

    return _is1, _is2
end

function peak_degree_lmn(evals, evecs, u, b; lmn=1, ub_pt=:up)
    degrees = zeros(Int,length(evals))
	if ub_pt == :up
		specid=1
	elseif ub_pt == :ut
		specid=2
	elseif ub_pt == :bp
		specid=3
	elseif ub_pt == :bt
		specid=4
	else
		error("specid must be :up, :ut, :bp, or :bt !")
	end
   spec = spectrum(evecs, u, b; lmn)[specid]
	@assert length(evals) == size(evecs,2)
    for i in axes(evecs,2)
        specs = @views spec[:,i]
        peak = findmax(specs)[2]
        degrees[i]=peak
    end
    return degrees
end

function energydiff(evecs, u, b, cutoff; thresh=1e-2, lmn=1, ub_pt=:up)
    ediff = zeros(Bool,size(evecs,2))
	if ub_pt == :up
		specid=1
	elseif ub_pt == :ut
		specid=2
	elseif ub_pt == :bp
		specid=3
	elseif ub_pt == :bt
		specid=4
	else
		error("specid must be :up, :ut, :bp, or :bt !")
	end
   spec = spectrum(evecs, u, b; lmn)[specid]
    for i in axes(evecs,2)
        spec_large = @views spec[1:(cutoff+1),i]
        spec_small = @views spec[(cutoff+2):end,i]
        if mean(spec_small)/mean(spec_large) < thresh
            ediff[i] = true
        end
    end
    return ediff
end

"""
$(TYPEDSIGNATURES)

Find observationally relevant solutions, based on the poloidal magnetic field.
"""
function observability_filter(evals, evecs, u, b;
	ωlow = 0.57, 
	ωhigh = 12.6, 
	Qlow = 1.0,
	lthresh = 17
	)
	
	@inline Q(f) = abs(imag(f)/2real(f))
	ω = imag.(evals)
	_ωfilter = @.(ωlow < abs(ω) < ωhigh)
	_Qfilter = Q.(evals) .> Qlow
	is = eachindex(evals)[_ωfilter .& _Qfilter]
	
	@views peakl = peak_degree_lmn(evals[is], evecs[:, is], u, b, lmn=1, ub_pt=:bp)
	# @views peakm = peak_degree_lmn(evals[is], evecs[:, is], u, b, lmn=2, ub_pt=:bp)
	@views peakn = peak_degree_lmn(evals[is], evecs[:, is], u, b, lmn=3, ub_pt=:bp)
	_lfilter = (peakl.<=lthresh) .& (peakn .<= lthresh÷2)
	
	is = is[_lfilter]
	return is
end

function observability_filter_ediff(evals, evecs, u, b;
	ωlow = 0.57, 
	ωhigh = 12.6, 
	Qlow = 1.0,
	lthresh = 17,
    thresh= 1e-2
	)
	
	@inline Q(f) = abs(imag(f)/2real(f))
	ω = imag.(evals)
	_ωfilter = @.(ωlow < abs(ω) < ωhigh)
	_Qfilter = Q.(evals) .> Qlow
	is = eachindex(evals)[_ωfilter .& _Qfilter]
	
	@views ediff_l = energydiff(evecs[:, is], u, b, lthresh; thresh, lmn=1, ub_pt=:bp)
	@views ediff_n = energydiff(evecs[:, is], u, b, lthresh÷2; thresh, lmn=3, ub_pt=:bp)
	# _lfilter = (peakl.<=lthresh) .& (peakn .<= lthresh÷2)
	
	is = is[ediff_l .& ediff_n]
	return is
end


end #module