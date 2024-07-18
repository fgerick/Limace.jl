module Processing

using Statistics
using ..Bases
using ..Bases: lmn_p, lmn_t



function getindices(u0,b0, u1,b1)
	lmnp0 = lmn_p(u0)
	lmnt0 = lmn_t(u0)
	np0 = length(lmnp0)
	nu = np0 + length(lmnt0)
	lmnpnew = lmn_p(u1)
	lmntnew = lmn_t(u1)
	indsp = findall(in(lmnpnew),lmnp0)
	indst = findall(in(lmntnew),lmnt0).+np0
	indsu = vcat(indsp,indst) 

	lmnp0 = lmn_p(b0)
	lmnt0 = lmn_t(b0)
	np0 = length(lmnp0)
	lmnpnew = lmn_p(b1)
	lmntnew = lmn_t(b1)
	indsp = findall(in(lmnpnew),lmnp0).+nu
	indst = findall(in(lmntnew),lmnt0).+np0.+nu
	indsb = vcat(indsp,indst)
	return vcat(indsu,indsb)
end


function usectionrev(::Val{true}, λ, λ2, evecs, evecs2, u0, b0, u1, b1;threshc=0.99, threshλ=0.1)
	@assert u1.N > u0.N
	nth = Threads.nthreads()
	is = [Int[] for _ in 1:nth]
	is2 = [Int[] for _ in 1:nth]
	inds = getindices(u0, b0, u1, b1)
	
	Threads.@threads for i in eachindex(λ)
		it = Threads.threadid()
		found = false
		ω = λ[i]
		for i2 in eachindex(λ2)
			ω2 = λ2[i2]
			if (!found) && (abs(ω)*(1-threshλ) < abs(ω2) < (1+threshλ)*abs(ω))
				c = cor(@views(evecs[:,i]), @views(evecs2[inds,i2]))
				if (abs(c)>threshc) && (sign(imag(ω)) == sign(imag(ω2)))
					push!(is[it], i)
					push!(is2[it],i2)
					found = true
					@goto λi
				end
			end
		end
		@label λi
	end
	return vcat(is...), vcat(is2...)
end

end #module