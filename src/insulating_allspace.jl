module InsulatingMFBasis

using SparseArrays

# The basis elements are normalized so that ∫₀^∞ B⋅B dV = 1.
# There is no external contribution for the toroidal components.
# For the poloidal components with n=1, there is an external contribution (r>1).
# To normalize the basis elements so that ∫₀¹ B⋅B dV = 1, the following norm function can be used:
function unitspherenorm(l,n)
    if n==1
        return sqrt((6+8l*(l+2))/(6+l*(6l+11)))
    else
        return 1.0
    end
end


#lmns

function lmn_bpol(N, ms = -N:N, ns = false) 
    if ns != false
        [(l,m,n) for l in 1:N for m in ms for n in ns if abs(m)<=l]
    else
        [(l,m,n) for l in 1:N for m in ms for n in 1:(N-l+1)÷2 if abs(m)<=l] 
        # [(l,m,n) for n in 1:N for l in 1:N for m in ms if (abs(m)<=l) && (n<=((N-l+1)÷2))] 
    end
end

function lmn_btor(N, ms = -N:N, ns = false) 
    if ns != false
        [(l,m,n) for l in 1:N  for m in ms for n in ns if abs(m)<=l]
    else
        [(l,m,n) for l in 1:N for m in ms for n in 1:((N-l)÷2)  if abs(m)<=l]  
        # [(l,m,n) for n in 1:N for l in 1:N for m in ms  if (abs(m)<=l) && (n<=((N-l)÷2))] 
    end
end

# function lmn_bpol(N, ms=-N:N, ns=0)
#     vcat(_lmn_bpol(1,ms,ns),[setdiff(_lmn_bpol(n,ms,ns),_lmn_bpol(n-1,ms,ns)) for n in 2:N]...)
# end
# function lmn_btor(N, ms=-N:N, ns=0)
#     vcat(_lmn_btor(1,ms,ns),[setdiff(_lmn_btor(n,ms,ns),_lmn_btor(n-1,ms,ns)) for n in 2:N]...)
# end
n(N) = ((-1)^(2*N)*(-1 + N)*N*(5 + 2*N))÷6
np(N) = ((-1)^N*(3 + (-1)^N*(-3 + 2*N*(-1 + N*(3 + N)))))÷12
nt(N) = ((-1)^N*(-3 + (-1)^N*(3 - 8*N + 2*N^3)))÷12

nlp(N,l) = ((-1)^N*(3 - 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(1 - 6*l - 4*l^2 + 6*(2 + l)*N))))÷12
nlt(N,l) = ((-1)^N*(-3 + 3*(-1)^l*(1 + l) + (-1)^N*(12*(-1 + (-1)^(2*l)) + l*(-11 - 4*l*(3 + l) + 12*N + 6*l*N))))÷12

lmn2k_p(l,m,n,N) = nlp(N,l-1) + (l+m)*((N-l+1)÷2) + n
lmn2k_t(l,m,n,N) = nlt(N,l-1) + (l+m)*((N-l)÷2) + n


function lmn_bpol_l(N, ms = -N:N, ns=0)
    lmn = lmn_bpol(N,ms,ns)
    lmnk = Vector{NTuple{4,Int}}[]
    L = N
    for _ in 1:L
        push!(lmnk,NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l,m,n = lmn[k]
        push!(lmnk[l], (k,l,m,n))
    end
    return lmnk
end

function lmn_btor_l(N, ms = -N:N, ns=0)
    lmn = lmn_btor(N,ms,ns)
    lmnk = Vector{NTuple{4,Int}}[]
    L = N
    for _ in 1:L
        push!(lmnk,NTuple{4,Int}[])
    end

    for k in eachindex(lmn)
        l,m,n = lmn[k]
        push!(lmnk[l], (k,l,m,n))
    end
    return lmnk
end


# function lmn_bpol_ml(N, ms = 0:(N-1), ns=0)
#     # lmn = [[(l,m,n) for m in ms for n in 1:((N-l+1)÷2+1) if abs(m)<=l] for l in 1:(N-1)] 
#     # if (ns != 0)
#     #     lmn = [[[(l,m,n) for n in ns if abs(m)<=l] for m in ms] for l in 1:(N-1)] 
#     # else
#     #     lmn = [[[(l,m,n) for n in 1:((N-l+1)÷2+1) if abs(m)<=l] for m in ms] for l in 1:(N-1)] 
#     # end
#     # lmnk = Vector{Vector{NTuple{4,Int}}}[]
#     # k=1
#     # for l in eachindex(lmn)
#     #     push!(lmnk,Vector{Vector{NTuple{4,Int}}}[])
#     #     for m in eachindex(lmn[l])
#     #         push!(lmnk[l],Vector{NTuple{4,Int}}[])
#     #         for (n,lmn) in enumerate(lmn[l][m])
#     #             push!(lmnk[l][m],(k,lmn...))
#     #             k+=1
#     #         end
#     #     end
#     # end

#     lmnk = Vector{Vector{NTuple{4,Int}}}[]
#     k=1
#     for l in 1:(N-1)
#         push!(lmnk,Vector{Vector{NTuple{4,Int}}}[])
#         for m in ms
#             if abs(m)<=l
#                 push!(lmnk[l],Vector{NTuple{4,Int}}[])
#                 for n in (ns==0 ? (1:((N-l+1)÷2 +1)) : ns)
#                     push!(lmnk[l][m+1],(k,l,m,n))
#                     k+=1
#                 end
#             end
#         end
#     end
#     return lmnk
# end

#inner products

@inline function _inner_tt!(is,js,aijs, i,j, l,n,n2)
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

@inline function _inner_ss!(is,js,aijs, i,j, l,n,n2)
    
   
    if n==n2 
        aij = one(l)
        #also add the all space contribution here (only nonzero for n==n2==1):
        # if n==1
        #     # aij += l^2*(1 + l)*(5 + 2l)^2
        #     aij += (l*(5 + 2*l))/(6 + l*(11 + 6*l))
        # end
        push!(is,i)
        push!(js,j)
        push!(aijs,aij)
    elseif (n==n2+1) 
        # if n2 == 1
        #     aij = -(((1 + 2l)*(9 + 2l))/(sqrt(2*(7 + 2l)*(9 + 2l)*(6 + l*(11 + 6l)))))
        # else
            aij = -sqrt(1 + 3/(5 - 2*l - 4*n) + 3/(-1 + 2*l + 4*n))/2
        # end
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


function lhs(N,m; ns = false, Ω::T = 1.0) where T
    lmn_p = lmn_bpol(N,m,ns)
    lmn_t = lmn_btor(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        _inner_ss!(is,js,aijs,i,i, T(l),T(n),T(n))
        if n>1
            _inner_ss!(is,js,aijs,i,i-1, T(l),T(n),T(n-1))
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        _inner_tt!(is,js,aijs, i+np,i+np, T(l),T(n),T(n))
        if n>1
            _inner_tt!(is,js,aijs, i+np,i-1+np, T(l),T(n),T(n-1))
        end
    end

    LHS = sparse(is,js,aijs, nu, nu)
    return LHS

end

#diffusion

@inline function _diffusion_tt!(is,js,aijs, i,j, l,m,n,n2; η = 1.0) 
    
   
    if n==n2 
        aij = -((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,η*aij)
    end
        
    return nothing
end


@inline function _diffusion_ss!(is,js,aijs, i,j, l,m,n,n2; η = 1.0) 
    
   
    if n==n2 
        # if n==1
        #     aij = -(5+7l+2l^2)/2
        # else
            aij = -((-3 + 2*l + 4*n)*(1 + 2*l + 4*n))/2
        # end
        push!(is,i)
        push!(js,j)
        push!(aijs,η*aij)
    end

    return nothing
end

function rhs_diffusion(N,m; ns = false, η::T = 1.0) where T
    lmn_p = lmn_bpol(N,m,ns)
    lmn_t = lmn_btor(N,m,ns)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = T.((l,m,n))
        _diffusion_ss!( is,js,aijs,i,i, l,m,n,n; η)
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = T.((l,m,n))
        _diffusion_tt!( is,js,aijs, i+np,i+np, l,m,n,n; η)  
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end



end