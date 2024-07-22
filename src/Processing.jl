module Processing

using SparseArrays
using Statistics
using ..Bases
using ..Bases: lmn_p, lmn_t, lmn_p_bc, lmn_t_bc, lmn2k_p_dict, lmn2k_t_dict



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

	lmnpnew = lmn_p_bc(u1)
	lmntnew = lmn_t_bc(u1)
	indsp = findall(in(lmnpnew),lmnp0)
	indst = findall(in(lmntnew),lmnt0).+np0
	indsu_bc = vcat(indsp,indst) 

	lmnp0 = lmn_p(b0)
	lmnt0 = lmn_t(b0)
	np0 = length(lmnp0)

	lmnpnew = lmn_p(b1)
	lmntnew = lmn_t(b1)
	indsp = findall(in(lmnpnew),lmnp0).+nu
	indst = findall(in(lmntnew),lmnt0).+np0.+nu
	indsb = vcat(indsp,indst)

	lmnpnew = lmn_p_bc(b1)
	lmntnew = lmn_t_bc(b1)
	indsp = findall(in(lmnpnew),lmnp0).+nu
	indst = findall(in(lmntnew),lmnt0).+np0.+nu
	indsb_bc = vcat(indsp,indst)



	return vcat(indsu_bc, indsb_bc), vcat(indsu,indsb)
end


function getindices_nobc(u,b)

	lmnpu = lmn_p_bc(u)
	lmnpud = lmn2k_p_dict(u)
	ipu = [lmnpud[i] for i in lmnpu]
	jpu = [lmnpud[i] for i in lmn_p(u)]
	np0 = length(lmn_p(u))

	lmntu = lmn_t_bc(u)
	lmntud = lmn2k_t_dict(u)
	itu = [lmntud[i]+np0 for i in lmntu]
	iu = vcat(ipu,itu)

	jtu = [lmntud[i]+np0 for i in lmn_t(u)]
	ju = vcat(jpu,jtu)

	nu = np0 + length(lmn_t(u))

	lmnpb = lmn_p_bc(b)
	lmnpbd = lmn2k_p_dict(b)
	ipb = [lmnpbd[i]+nu for i in lmnpb]
	jpb = [lmnpbd[i]+nu for i in lmn_p(b)]
	lmntb = lmn_t_bc(b)
	lmntbd = lmn2k_t_dict(b)

	np0 = length(lmn_p(b))
	itb = [lmntbd[i]+np0+nu for i in lmntb]
	jtb = [lmntbd[i]+np0+nu for i in lmn_t(b)]
	ib = vcat(ipb,itb)
	jb = vcat(jpb,jtb)

	return vcat(iu,ib), vcat(ju,jb)
end

function reducematrix(a,u0,b0,u1,b1)
	i,j=getindices(u0,b0, u1,b1)	
	n = length(u1)+length(b1)
	out = spzeros(ComplexF64,n,n)
	a1 = a[i,j]
	i2,j2 = getindices_nobc(u1,b1)
	out[i2,j2] = a1
	return out
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