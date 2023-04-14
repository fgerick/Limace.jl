using Limace.DiscretePart: s_mf_pre, d_s_mf_pre, d2_s_mf_pre, d3_s_mf_pre
using Limace.DiscretePart: t_mf_pre, d_t_mf_pre, d2_t_mf_pre, d3_t_mf_pre
using Limace.DiscretePart: s_in_pre, d_s_in_pre, d2_s_in_pre, d3_s_in_pre
using Limace.DiscretePart: t_in_pre, d_t_in_pre, d2_t_in_pre, d3_t_in_pre
using Limace.DiscretePart: t_in, s_in, t_mf, s_mf, rquad, jacobis
const DP = Limace.DiscretePart


@testset "lorentz" begin
	N = 50
	r,wr = rquad(N)
	lmna = (1,1,1)
	lmns = [(l,m,n) for l in 1:3 for m in -l:l for n in [1]]
	lmnb = (1,1,1)

	js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
	js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]
	rls = [r.^l for l in 1:N]

	Base.@propagate_inbounds Smf(l,m,n,r,i) = s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds dSmf(l,m,n,r,i) = d_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds d2Smf(l,m,n,r,i) = d2_s_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds d3Smf(l,m,n,r,i) = d3_s_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds Tmf(l,m,n,r,i) = t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds dTmf(l,m,n,r,i) = d_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds d2Tmf(l,m,n,r,i) = d2_t_mf_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds d3Tmf(l,m,n,r,i) = d3_t_mf_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds tu(l,m,n,r,i) = t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds dtu(l,m,n,r,i) = d_t_in_pre(js_a0,rls,l,m,n,r,i)
	Base.@propagate_inbounds d2tu(l,m,n,r,i) = d2_t_in_pre(js_a0,rls,l,m,n,r,i)

	Base.@propagate_inbounds su(l,m,n,r,i) = s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds dsu(l,m,n,r,i) = d_s_in_pre(js_a1,rls,l,m,n,r,i)
	Base.@propagate_inbounds d2su(l,m,n,r,i) = d2_s_in_pre(js_a1,rls,l,m,n,r,i)

	for lmna in lmns, lmnc in lmns
		l1 = DP._induction_sSS_pre(lmna, lmnb, lmnc, r,wr, su,Smf,Smf,  dsu, dSmf, dSmf, d2su, d2Smf)
		l2 = DP._induction_sSS(lmna, lmnb, lmnc, r,wr, s_in, s_mf, s_mf)
		@test l1 ≈ l2 atol=1e-14

		l1 = DP._induction_sTS_pre(lmna, lmnb, lmnc, r,wr, su,Tmf,Smf,dsu, dTmf, dSmf)
		l2 = DP._induction_sTS(lmna, lmnb, lmnc, r,wr, s_in,t_mf,s_mf)
		@test l1 ≈ l2 atol=1e-14

		l1 = DP._induction_tSS_pre(lmna, lmnb, lmnc, r,wr, tu, Smf,Smf, dtu, dSmf, dSmf, t_in, s_mf, s_mf)
		l2 = DP._induction_tSS(lmna, lmnb, lmnc, r,wr, t_in,s_mf,s_mf)
		@test l1 ≈ l2 atol=1e-14
		
	end
	# for lmna in lmns, lmnc in lmns
	# 	l1 = DP._lorentz_SSs_pre(lmna, lmnb, lmnc, r,wr, Smf,Smf,su, dSmf, d2Smf, d3Smf, dSmf )
	# 	l2 = DP._lorentz_SSs(lmna, lmnb, lmnc, r,wr, s_mf, s_mf, s_in)
	# 	@test l1 ≈ l2 atol=1e-14
		
	# 	l1 = DP._lorentz_STs_pre(lmna, lmnb, lmnc, r,wr, Smf,Tmf,su, dSmf, d2Smf, dTmf, d2Tmf)
	# 	l2 = DP._lorentz_STs(lmna, lmnb, lmnc, r,wr, s_mf, t_mf, s_in)
	# 	@test l1 ≈ l2 atol=1e-14

	# 	l1 = DP._lorentz_TTs_pre(lmna, lmnb, lmnc, r,wr, Tmf,Tmf,su, dTmf, dTmf)	
	# 	l2 = DP._lorentz_TTs(lmna, lmnb, lmnc, r,wr, t_mf,t_mf,s_in)	
	# 	@test l1 ≈ l2 atol=1e-14

	# 	l1 = DP._lorentz_STt_pre(lmna, lmnb, lmnc, r,wr, Smf,Tmf,tu, dSmf, dTmf)
	# 	l2 = DP._lorentz_STt(lmna, lmnb, lmnc, r,wr, s_mf,t_mf,t_in)
	# 	@test l1 ≈ l2 atol=1e-14

	# 	l1 = DP._lorentz_TTt_pre(lmna, lmnb, lmnc, r,wr, Tmf, Tmf, tu)
	# 	l2 = DP._lorentz_TTt(lmna, lmnb, lmnc, r,wr, t_mf, t_mf, t_in)
	# 	@test l1 ≈ l2 atol=1e-14
	# end
	# N = 10
	# r,wr = rquad(2N)
	# js_a1 = [jacobis(N,1.0,l+1/2,r) for l in 1:N]
	# js_a0 = [jacobis(N,0.0,l+1/2,r) for l in 1:N]

	# @test DP.rhs_lorentz_bpol_pre(N,-N:N, lmnb, r, wr, js_a1, js_a0) ≈ DP.rhs_lorentz_bpol(N,-N:N, lmnb)
end