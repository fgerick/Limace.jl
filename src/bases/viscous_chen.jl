module ChenBasis

using SparseArrays

function lmn_upol(N, ms = 0:N, ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:(N-1) for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:(N-1) for n in 1:((N-l+1)÷2+1) if abs(m)<=l] 
    end
end

function lmn_utor(N, ms = 0:N, ns = 0) 
    if (ns != 0)
        [(l,m,n) for m in ms for l in 1:N for n in ns if abs(m)<=l]
    else
        [(l,m,n) for m in ms for l in 1:N for n in 1:((N-l)÷2+1) if abs(m)<=l] 
    end
end

function lmn_upol_l(N, ms = 0:N, ns=0)
    # lmn = [[(l,m,n) for m in ms for n in 1:((N-l+1)÷2+1) if abs(m)<=l] for l in 1:(N-1)] 
    if (ns != 0)
        lmn = [[(l,m,n) for m in ms for n in ns if abs(m)<=l] for l in 1:(N-1)]
    else
        lmn = [[(l,m,n) for m in ms for n in 1:((N-l+1)÷2+1) if abs(m)<=l] for l in 1:(N-1)]
    end
    lmnk = Vector{NTuple{4,Int}}[]
    k=1
    for l in eachindex(lmn)
        push!(lmnk,NTuple{4,Int}[])
        for lmn in lmn[l]
            push!(lmnk[l], (k,lmn...))
            k+=1
        end
    end
    return lmnk
end

function lmn_utor_l(N, ms = 0:N, ns=0)
    if (ns != 0)
        lmn = [[(l,m,n) for m in ms for n in ns if abs(m)<=l] for l in 1:N] 
    else
        lmn = [[(l,m,n) for m in ms for n in 1:((N-l)÷2+1) if abs(m)<=l]  for l in 1:N]
    end
    lmnk = Vector{NTuple{4,Int}}[]
    k=1
    for l in eachindex(lmn)
        push!(lmnk,NTuple{4,Int}[])
        for lmn in lmn[l]
            push!(lmnk[l], (k,lmn...))
            k+=1
        end
    end
    return lmnk
end

@inline function _inner_tt(is,js,aijs, i,j, l,n,n2)
    
   
    if n==n2 
        aij = one(l)
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    elseif (n==n2+1)
        # n = n-1
        aij = -sqrt(1 - 3/(1 + 2*l + 4*(-1 + n)) + 3/(5 + 2*l + 4*(-1 + n)))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
        #symmetric
        push!(is,j)
        push!(js,i)
        push!(aijs,aij)

    end
        
    return nothing
end

@inline function _inner_ss(is,js,aijs, i,j, l,n,n2)
    
   
    if n==n2 
        aij = one(l)
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    elseif (n==n2+1) 
        # n = n-1
        aij = -1/2sqrt(1 + 3/(2 + 4*l + 8*(-1 + n)) - 3/(18 + 4*l + 8*(-1 + n)))
        # aij = -1/(2*sqrt(1 + 3/(2 + 4*l + 8*n) - 3/(18 + 4*l + 8*n)))
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
        #symmetric
        push!(is,j)
        push!(js,i)
        push!(aijs,aij)
    end
        
    return nothing
end

# function lhs_original(N,m; Ω::T = 1.0) where T
#     lmn_p = lmn_upol(N,m)
#     lmn_t = lmn_utor(N,m)

#     np = length(lmn_p)
#     @show np
#     nt = length(lmn_t)
#     nu = np+nt

#     is,js,aijs = Int[], Int[], Complex{T}[]


#     for (i,(l,m,n)) in enumerate(lmn_p)
#         # l,m,n = T.((l,m,n))
#         for (j, (l2,m2,n2)) in enumerate(lmn_p)
#             # l2,m2,n2 = T.((l2,m2,n2))
#             if (l==l2) && (m==m2)
#                 _inner_ss(is,js,aijs,i,j, T(l),T(n),T(n2))
#             # j+=1
#             end
#         end
#     end

#     for (i,(l,m,n)) in enumerate(lmn_t)
#         for (j, (l2,m2,n2)) in enumerate(lmn_t)
#             if (l==l2) && (m==m2)
#                 # l2,m2,n2 = T.((l2,m2,n2))
#                 _inner_tt(is,js,aijs, i+np,j+np, T(l),T(n),T(n2))
#                 # j+=1
#             end
#         end
#     end

#     LHS = sparse(is,js,aijs, nu, nu)
#     return LHS

# end


function lhs(N,m; ns = 0, Ω::T = 1.0) where T
    lmn_p = lmn_upol(N,m,ns)
    lmn_t = lmn_utor(N,m,ns)

    np = length(lmn_p)
    # @show np
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        _inner_ss(is,js,aijs,i,i, T(l),T(n),T(n))
        if n>1
            _inner_ss(is,js,aijs,i,i-1, T(l),T(n),T(n-1))
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        _inner_tt(is,js,aijs, i+np,i+np, T(l),T(n),T(n))
        if n>1
            _inner_tt(is,js,aijs, i+np,i-1+np, T(l),T(n),T(n-1))
        end
    end

    LHS = sparse(is,js,aijs, nu, nu)
    return LHS

end


@inline function _coriolis_tt(is,js,aijs, i,j, l,m,n,n2; Ω = 2) 
    if n==n2
        aij = Ω*im*m/(2l*(l+1))
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    elseif (n==n2+1)
        n = n-1
        aij = -im*m/(2l*(1 + l)*sqrt(1 + 3/(2*(-1 + 2*l + 4*n)) - 3/(14 + 4*l + 8*n)))
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    elseif (n==n2-1)
        n = n+1
        aij = -im*m*sqrt(1 - 3/(-3 + 2*l + 4*n) + 3/(1 + 2*l + 4*n))/(2l*(1 + l))
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    end

end
@inline function _coriolis_ss(is,js,aijs, i,j, l,m,n,n2; Ω = 2) 
    if n==n2
        aij = Ω*im*m/(2l*(l+1))
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    elseif (n==n2+1)
        n = n-1
        aij = -im*m/(2l*(1 + l)*sqrt(1 + 3/(2 + 4*l + 8*n) - 3/(18 + 4*l + 8*n)))
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    elseif (n==n2-1)
        n=n+1
        aij = -im*m*sqrt(1 - 3/(-1 + 2*l + 4*n) + 3/(3 + 2*l + 4*n))/(2l*(1 + l))    
        push!(is,i)
        push!(js,j)
        push!(aijs,Ω*aij)
    end

end


@inline function _coriolis_ts(is,js,aijs, i,j, l,l2,m,m2,n,n2; Ω = 2) 
    if (m != m2)
        return nothing
    end
    if (l==l2+1)
        if n==n2
            aij = ((-1 + l^2)*sqrt(((l - m)*(l + m))/(1 - 5*l^2 + 4*l^4)))/l
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif (n == n2-1)
            aij = -((-1 + l^2)*sqrt(((l - m)*(l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((1 - 5*l^2 + 4*l^4)*(1 + 2*l + 4*n)*(5 + 2*l + 4*n))))/l/2 
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif (n == n2+1)
            aij = -((-1 + l^2)*sqrt((l - m)*(l + m))*sqrt(-5 + 2*l + 4*n))/(l*sqrt(((1 - 5*l^2 + 4*l^4)*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        end
    elseif (l==l2-1)
        if n==n2
            aij = -sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((3 + 4*l*(2 + l))*(1 + 2*l + 4*n)*(5 + 2*l + 4*n)))/(1 + l)/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif n==n2+1
            aij = sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m))/(3 + 4*l*(2 + l)))/(1 + l)
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif n==n2+2
            aij = -(sqrt(l*(2 + l)*(1 + l - m)*(1 + l + m))*sqrt(-5 + 2*l + 4*n))/((1 + l)*sqrt(((3 + 4*l*(2 + l))*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        end
       
    end
    return nothing
end

@inline function _coriolis_st(is,js,aijs, i,j, l2,l,m2,m,n2,n; Ω = 2) 
    if (m != m2)
        return nothing
    end
    if (l==l2+1)
        if n==n2
            aij = -((-1 + l^2)*sqrt(((l - m)*(l + m))/(1 - 5*l^2 + 4*l^4)))/l
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif (n == n2-1)
            aij = ((-1 + l^2)*sqrt(((l - m)*(l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((1 - 5*l^2 + 4*l^4)*(1 + 2*l + 4*n)*(5 + 2*l + 4*n))))/l/2 
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif (n == n2+1)
            aij = ((-1 + l^2)*sqrt((l - m)*(l + m))*sqrt(-5 + 2*l + 4*n))/(l*sqrt(((1 - 5*l^2 + 4*l^4)*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        end
    elseif (l==l2-1)
        if n==n2
            aij = sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m)*(-1 + 2*l + 4*n)*(7 + 2*l + 4*n))/((3 + 4*l*(2 + l))*(1 + 2*l + 4*n)*(5 + 2*l + 4*n)))/(1 + l)/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif n==n2+1
            aij = -sqrt((l*(2 + l)*(1 + l - m)*(1 + l + m))/(3 + 4*l*(2 + l)))/(1 + l)
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        elseif n==n2+2
            aij = (sqrt(l*(2 + l)*(1 + l - m)*(1 + l + m))*sqrt(-5 + 2*l + 4*n))/((1 + l)*sqrt(((3 + 4*l*(2 + l))*(-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/(3 + 2*l + 4*n)))/2
            push!(is,i)
            push!(js,j)
            push!(aijs,Ω*aij)
        end
       
    end
    return nothing
end


@inline function _viscous_tt(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
    if n==n2 
        aij = -((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,ν*aij)
    end
        
    return nothing
end


@inline function _viscous_ss(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
    if n==n2 
        aij = -((1 + 2*l + 4*n)*(5 + 2*l + 4*n))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,ν*aij)
    end
    # end
        
    return nothing
end

function rhs(N,m; Ω::T = 2.0, ν::T = 1.0) where T
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    @show np
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = T.((l,m,n))
        for (j, (l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.((l2,m2,n2))
            if (l==l2) && (m==m2)
                _coriolis_ss(is,js,aijs,i,j, l,m,n,n2; Ω)
                _viscous_ss( is,js,aijs,i,j, l,m,n,n2; ν)
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_st(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = T.((l,m,n))
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
                l2,m2,n2 = T.((l2,m2,n2))
                _coriolis_tt(is,js,aijs, i+np,j+np, l,m,n,n2; Ω)
                _viscous_tt( is,js,aijs, i+np,j+np, l,m,n,n2; ν)  
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_ts(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

@inline function rhs_coriolis_1(indices,lmn_p, lmn_t; Ω::T = 2.0 ) where T

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for i in indices
        # (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = T.(lmn_p[i])
        for (j, (l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.((l2,m2,n2))
            if (l==l2) && (m==m2)
                _coriolis_ss(is,js,aijs,i,j, l,m,n,n2; Ω)
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_st(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end

    return is, js, aijs
end


@inline function rhs_coriolis_2(indices,lmn_p, lmn_t; Ω::T = 2.0 ) where T

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for i in indices
        l,m,n = T.(lmn_t[i])
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
            l2,m2,n2 = T.(lmn_t[j])
            _coriolis_tt(is,js,aijs, i+np,j+np, l,m,n,n2; Ω)
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.(lmn_p[j])
            _coriolis_ts(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end

    return is, js, aijs
end


@inline function _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0 ) where T

    is,js,aijs = Int[], Int[], Complex{T}[]

    for (l,ilmn) in enumerate(lmn_p_l)
        for (i,_,m,n) in ilmn
            for (j,l2,m2,n2) in lmn_p_l[l]
                if (m==m2)
                    _coriolis_ss(is,js,aijs, i,j, T(l),T(m),T(n),T(n2); Ω)
                end
            end
            lrange = (l== 1) ? [2] : ((l==N) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_t_l,lrange)
                # l2,m2,n2 = T.(lmn_p[j])
                for (j,l2,m2,n2) in lmn2
                    _coriolis_st(is,js,aijs, i,j+np, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                end
            end
        end
    end

    return is, js, aijs
end
@inline function _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω::T = 2.0 ) where T

    # np = length(lmn_p)
    # nt = length(lmn_t)
    # nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (l,ilmn) in enumerate(lmn_t_l)
        for (i,_,m,n) in ilmn
            # l,m,n = T.((l,m,n))
            for (j,l2,m2,n2) in lmn_t_l[l]
                if (m==m2)
                # l2,m2,n2 = T.(lmn_t[j])
                _coriolis_tt(is,js,aijs, i+np,j+np, T(l),T(m),T(n),T(n2); Ω)
                end
            end
            # lrange = l == 1 ? [1,2] : (l-1:l+1)
            lrange = (l== 1) ? [2] : ((l>=N-1) ? [l-1] : [l-1,l+1])
            for lmn2 in view(lmn_p_l,lrange)
                # l2,m2,n2 = T.(lmn_p[j])
                for (j,l2,m2,n2) in lmn2
                    _coriolis_ts(is,js,aijs, i+np,j, T(l),T(l2),T(m),T(m2),T(n),T(n2); Ω)
                end
            end
        end
    end

    return is, js, aijs
end
function rhs_coriolis(N,m; ns = 0, Ω::T = 2.0) where T
    lmn_p = lmn_upol(N,m,ns)
    lmn_t = lmn_utor(N,m,ns)
    lmn_p_l = lmn_upol_l(N,m,ns)
    lmn_t_l = lmn_utor_l(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    # is, js, aijs = rhs_coriolis_1(eachindex(lmn_p), lmn_p, lmn_t; Ω )
    # is2, js2, aijs2 = rhs_coriolis_2(eachindex(lmn_t), lmn_p, lmn_t; Ω )
    is, js, aijs = _rhs_coriolis_1(N, np, lmn_p_l, lmn_t_l; Ω )
    is2, js2, aijs2 = _rhs_coriolis_2(N, np, lmn_p_l, lmn_t_l; Ω )

    append!(is,is2)
    append!(js,js2)
    append!(aijs,aijs2)

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end


function rhs_viscosity(N,m; ns = 0, ν::T = 1.0) where T
    lmn_p = lmn_upol(N,m,ns)
    lmn_t = lmn_utor(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = T.((l,m,n))
        # for (j, (l2,m2,n2)) in enumerate(lmn_p)
        #     l2,m2,n2 = T.((l2,m2,n2))
        #     if (l==l2) && (m==m2) && (n==n2)
        _viscous_ss( is,js,aijs,i,i, l,m,n,n; ν)
        #     end
        # end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = T.((l,m,n))
        # for (j, (l2,m2,n2)) in enumerate(lmn_t)
        #     if (l==l2) && (m==m2) && (n==n2)
        #         l2,m2,n2 = T.((l2,m2,n2))
        _viscous_tt( is,js,aijs, i+np,i+np, l,m,n,n; ν)  
        #     end
        # end
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

end