function _inertial_SS(lmna, lmnb, r,wr, Sa,Sb)
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _Sa = r->Sa(la,ma,na,r)
    @inline _Sb = r->Sb(lb,mb,nb,r)

    @inline f = r-> inners(_Sa,_Sb, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end

function _inertial_TT(lmna, lmnb, r,wr, Ta,Tb)
    la,ma,na = lmna
    lb,mb,nb = lmnb

    @inline _Ta = r->Ta(la,ma,na,r)
    @inline _Tb = r->Tb(lb,mb,nb,r)

    @inline f = r-> innert(_Tb,_Ta, lb, r)

    aij = ∫dr(f,r,wr)
    return aij
end


function lhs_inertial_b(N,m; ns = false, η::T=1.0, thresh = sqrt(eps()), 
                            smf = s_mf, 
                            tmf = t_mf,
                            lmn_bp = InsulatingMFBasis.lmn_bpol(N,m,ns),
                            lmn_bt = InsulatingMFBasis.lmn_btor(N,m,ns),
                        ) where T
    np = length(lmn_bp)

    is,js,aijs = Int[],Int[],Complex{T}[]

    # r, wr = rquad(N+lu0+nu0+5)
    rwrs = [rquad(n+1) for n in 1:N]

    for (i,lmni) in enumerate(lmn_bp)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bp)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
			if (li==lj) && (mi==mj)
				r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
				aij = _inertial_SS(lmni,lmnj, r, wr, smf, smf)
				appendit!(is,js,aijs, i,j, aij; thresh)
			end
        end
    end

    for (i,lmni) in enumerate(lmn_bt)
        li,mi,ni = lmni
        for (j, lmnj) in enumerate(lmn_bt)
            lj,mj,nj = lmnj
            # !ncondition(lu0,ni,nu0,nj) && continue
			if (li==lj) && (mi==mj)
				r,wr = rwrs[min(N,li÷2+ni+lj÷2+nj+1)]
				aij = _inertial_TT(lmni, lmnj, r, wr, tmf, tmf)
				appendit!(is,js,aijs, i+np,j+np, aij; thresh)
			end
        end
    end
    nmatb = length(lmn_bp)+length(lmn_bt)

    return sparse(is,js,aijs,nmatb, nmatb)
end