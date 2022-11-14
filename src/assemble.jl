function lhs(N,m; thresh=1000eps())

    lns = [(l,n) for l=1:N for n in 1:(N-l+1)÷2]
    _N = length(lns) 
    # LHS = spzeros(_N,_N)
    is,js,aijs = Int[],Int[],Float64[]
    r, wr = gausslegendre(N+5)

    r = @. (r+1)/2
    wr/=2

    for (i,(l,n)) in enumerate(lns)
        # aij_s = dot(@.((s_visc(l,m,n_i,r)*s_visc(l,m,n_j,r)/r^2)*r^2),wr)
        aij_t = 0.0
        aij_s = 0.0
        for (r,wr) in zip(r,wr)
            aij_t += t_visc(l,m,n,r)^2*r^2*wr
            aij_s += (s_visc(l,m,n,r)^2*l*(l+1) + drs_visc(l,m,n,r)^2)*wr
        end
        aij_t*=l*(l+1)
        aij_s*=l*(l+1)
        if abs(aij_s)>thresh
            push!(is,i)
            push!(js,i)
            push!(aijs,aij_s)
        end
        if abs(aij_t)>thresh
            push!(is,i+_N)
            push!(js,i+_N)
            push!(aijs,aij_t)
        end
    end

    LHS = sparse(is,js,aijs)
    return LHS
end

C(l,m) = (l^2-1)*sqrt((l^2-m^2)/(4l^2-1))

function rhs(N,m; thresh=1000eps())

    lns = [(l,n) for l=1:N for n in 1:(N-l+1)÷2]
    _N = length(lns) 
    # LHS = spzeros(_N,_N)
    is,js,aijs = Int[],Int[],ComplexF64[]
    r, wr = gausslegendre(N+5)

    r = @. (r+1)/2
    wr/=2

    for (i,(l_i,n_i)) in enumerate(lns), (j,(l_j,n_j)) in enumerate(lns)
        if l_i == l_j
            l = l_i


            aij_tt = 0.0
            aij_ss = 0.0
            if n_i == n_j
                n = n_i
                for (r,wr) in zip(r,wr)
                    aij_tt += im*m/(l*(l+1))*t_visc(l,m,n,r)^2*r^2*wr
                    aij_ss += m/(l*(l+1))*(s_visc(l,m,n,r)^2*l*(l+1) + drs_visc(l,m,n,r)^2)*wr
                end
            end
            if abs(aij_ss)>thresh
                push!(is,i)
                push!(js,j)
                push!(aijs,aij_ss)
            end

            if abs(aij_tt)>thresh
                push!(is,i+_N)
                push!(js,j+_N)
                push!(aijs,aij_tt)
            end

        elseif l_i == l_j - 1
            l = l_i

            aij_ts = 0.0
            aij_st = 0.0
            if true
                n = n_i
                for (r,wr) in zip(r,wr)
                    #d_l^(l-1) = d_r + (l+1)/r
                    aij_ts += t_visc(l,m,n,r)/(l*(l+1))*r^2*wr * C(l,m)* (ds_visc(l-1,m,n,r) + (l+1)/r*s_visc(l-1,m,n,r))
                    aij_st -= s_visc(l,m,n,r)/(l*(l+1))*r^2*wr * C(l,m)* (dt_visc(l-1,m,n,r) + (l+1)/r*t_visc(l-1,m,n,r))
                end
            end
            if abs(aij_st)>thresh
                push!(is,i)
                push!(js,j+_N)
                push!(aijs,aij_st)
            end

            if abs(aij_ts)>thresh
                push!(is,i+_N)
                push!(js,j)
                push!(aijs,aij_ts)
            end
        elseif l_i == l_j + 1
            l = l_i

            aij_ts = 0.0
            aij_st = 0.0
            if true
                n = n_i
                for (r,wr) in zip(r,wr)
                    #d_l^(l+1) = d_r -l/r
                    aij_ts += t_visc(l,m,n,r)*r^2*wr * C(l+1,m)* (ds_visc(l+1,m,n,r) - l/r*s_visc(l+1,m,n,r))
                    aij_st -= s_visc(l,m,n,r)*r^2*wr * C(l+1,m)* (dt_visc(l+1,m,n,r) - l/r*t_visc(l+1,m,n,r))
                end
            end
            if abs(aij_st)>thresh
                push!(is,i)
                push!(js,j+_N)
                push!(aijs,aij_st)
            end

            if abs(aij_ts)>thresh
                push!(is,i+_N)
                push!(js,j)
                push!(aijs,aij_ts)
            end
        end
    end

    LHS = sparse(is,js,aijs)
    return LHS
end

#d_n^(n-1) = d_r + (n+1)/r
#d_n^(n+1) = d_r -n/r