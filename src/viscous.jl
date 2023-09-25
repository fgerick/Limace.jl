module ViscousBasis


using SparseArrays
using Distributed
using DistributedArrays
using Random

lmn_upol(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:(N-1) for n in 0:(N-l+1)÷2 if abs(m)<=l]
# lmn_upol(N, ms = 0:N) = [(l,m,n) for m in ms for n in 0:(N-1) for l in 1:(N-1) if (abs(m)<=l) && (n<=(N-l+1)÷2)]
lmn_utor(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:N for n in 0:((N-l)÷2) if abs(m)<=l]

# lmn_utor(N, ms = 0:N) = [(l,m,n) for m in ms for n in 0:N for l in 1:N if (abs(m)<=l) && (n<=(N-l)÷2)]

_coriolis_tt(l,m; Ω = 2) = Ω*im*m/(l*(l+1))
_coriolis_ss(l,m; Ω = 2) = _coriolis_tt(l,m; Ω)

# function _coriolis_ts(is,js,aijs, i,j, l,l2,m,m2,n,n2; Ω = 2) 
#     if (m != m2)
#         return nothing
#     end
#     if (l==l2+1) && (n>=n2-1)
#         if n+1==n2
#             aij = sqrt(l^2-1)*sqrt((l-m)*(l+m)*(1+n)*(3+2l+2n))/l/sqrt((4l^2-1)*(2+n)*(5+2l+2n))
#         else
#             aij = -sqrt((l^2-1)*(l-m)*(l+m)*(3+2l+4n2)*(7+2l+4n))/sqrt((4l^2-1)*(n+1)*(n+2)*(3+2l+2n)*(5+2l+2n))/l
#         end
#         push!(is,i)
#         push!(js,j)
#         push!(aijs,Ω*aij)
#     elseif (l==l2-1) && (n>=n2)
#         if n==n2
#             aij = (l+2)/(l+1)*sqrt(l*(n+1)*(l-m+1)*(l+m+1)*(2l+2n+3))/sqrt((l+2)*(2l+1)*(2l+3)*(n+2)*(2l+2n+5))
#         else
#             aij = -sqrt((2l+4n2+7)*l*(l+2)*(l-m+1)*(l+m+1)*(2l+4n+7))/((l+1)*sqrt((4l*(l+2)+3)*(n+1)*(n+2)*(2l+2n+3)*(2l+2n+5)))
#         end
#         push!(is,i)
#         push!(js,j)
#         push!(aijs,Ω*aij)
#     end
#     return nothing
# end

# function _coriolis_st(is,js,aijs, i,j, l2,l,m2,m,n2,n; Ω = 2) 
#     if (m != m2)
#         return nothing
#     end
#     if (l==l2+1) && (n>=n2-1)
#         if n+1==n2
#             aij = -sqrt(l^2-1)*sqrt((l-m)*(l+m)*(1+n)*(3+2l+2n))/l/sqrt((4l^2-1)*(2+n)*(5+2l+2n))
#         else
#             aij = sqrt((l^2-1)*(l-m)*(l+m)*(3+2l+4n2)*(7+2l+4n))/sqrt((4l^2-1)*(n+1)*(n+2)*(3+2l+2n)*(5+2l+2n))/l
#         end
#         push!(is,i)
#         push!(js,j)
#         push!(aijs,Ω*aij)
#     elseif (l==l2-1) && (n>=n2)
#         if n==n2
#             aij = -(l+2)/(l+1)*sqrt(l*(n+1)*(l-m+1)*(l+m+1)*(2l+2n+3))/sqrt((l+2)*(2l+1)*(2l+3)*(n+2)*(2l+2n+5))
#         else
#             aij = sqrt((2l+4n2+7)*l*(l+2)*(l-m+1)*(l+m+1)*(2l+4n+7))/((l+1)*sqrt((4l*(l+2)+3)*(n+1)*(n+2)*(2l+2n+3)*(2l+2n+5)))
#         end
#         push!(is,i)
#         push!(js,j)
#         push!(aijs,Ω*aij)
#     end
#     return nothing
# end


function _viscous_tt(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
    if n>=n2 
        aij = -((1 + n2)*sqrt(2 + n2)*sqrt((7 + 2*l + 4*n)*(3 + 2*l + 2*n2))*sqrt(5 + 2*l + 2*n2)*sqrt(7 + 2*l + 4*n2)*(9 + 2*n2*(7 + 2*n2) + l*(6 + 4*n2)))/sqrt((1 + n)*(2 + n)*(3 + 2*l + 2*n)*(5 + 2*l + 2*n)*(1 + n2))/6
    else
        #symmetric.
        nt = n
        n = n2
        n2 = nt
        aij = -((1 + n2)*sqrt(2 + n2)*sqrt((7 + 2*l + 4*n)*(3 + 2*l + 2*n2))*sqrt(5 + 2*l + 2*n2)*sqrt(7 + 2*l + 4*n2)*(9 + 2*n2*(7 + 2*n2) + l*(6 + 4*n2)))/sqrt((1 + n)*(2 + n)*(3 + 2*l + 2*n)*(5 + 2*l + 2*n)*(1 + n2))/6
    end
    push!(is,i)
    push!(js,j)
    push!(aijs,ν*aij)
    # end
        
    return nothing
end


# function _viscous_ss(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
#     if n<n2 
#         aij = -(n - n2)*(1 + n2)*(5 + 2*l + 2*n + 2*n2)*sqrt(((5 + 2*l + 4*n)*(5 + 2*l + 4*n2))/(1 + n2)^2)
#         push!(is,i)
#         push!(js,j)
#         push!(aijs,ν*aij)
#     end
#     # end
        
#     return nothing
# end

function rhs1(indices, lmn_p,lmn_t,Ω::T,ν::T) where T
    is,js,aijs = Int[], Int[], Complex{T}[]
    np = length(lmn_p)

    for i in indices
        l,m,n = T.(lmn_p[i])
        push!(is,i)
        push!(js,i)
        push!(aijs,_coriolis_ss(l,m; Ω))
        for (j, (l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.((l2,m2,n2))
            if (l==l2) && (m==m2)
                _viscous_ss(is,js,aijs,i,j, l,m,n,n2; ν)
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_st(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end
    return is, js, aijs
end

function rhs2(indices,lmn_p,lmn_t,Ω::T,ν::T) where T
    is,js,aijs = Int[], Int[], Complex{T}[]
    np = length(lmn_p)

    for i in indices
        l,m,n = T.(lmn_t[i])
        push!(is,i+np)
        push!(js,i+np)
        push!(aijs,_coriolis_tt(l,m; Ω))
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
                l2,m2,n2 = T.((l2,m2,n2))
                _viscous_tt(is,js,aijs, i+np,j+np, l,m,n,n2; ν)  
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = T.((l2,m2,n2))
            _coriolis_ts(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end
    
    return is, js, aijs
end
function rhs(N,m; Ω::T = 2.0, ν::T = 1.0) where T
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    # @show np
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], Complex{T}[]


    is, js, aijs = rhs1(eachindex(lmn_p), lmn_p, lmn_t, Ω, ν)
    # for (i,(l,m,n)) in enumerate(lmn_p)
    #     l,m,n = T.((l,m,n))
    #     push!(is,i)
    #     push!(js,i)
    #     push!(aijs,_coriolis_ss(l,m; Ω))
    #     for (j, (l2,m2,n2)) in enumerate(lmn_p)
    #         l2,m2,n2 = T.((l2,m2,n2))
    #         if (l==l2) && (m==m2)
    #             _viscous_ss(is,js,aijs,i,j, l,m,n,n2; ν)
    #         end
    #     end
    #     for (j,(l2,m2,n2)) in enumerate(lmn_t)
    #         l2,m2,n2 = T.((l2,m2,n2))
    #         _coriolis_st(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
    #     end
    # end

    # for (i,(l,m,n)) in enumerate(lmn_t)
    #     l,m,n = T.((l,m,n))
    #     push!(is,i+np)
    #     push!(js,i+np)
    #     push!(aijs,_coriolis_tt(l,m; Ω))
    #     for (j, (l2,m2,n2)) in enumerate(lmn_t)
    #         if (l==l2) && (m==m2)
    #             l2,m2,n2 = T.((l2,m2,n2))
    #             _viscous_tt(is,js,aijs, i+np,j+np, l,m,n,n2; ν)  
    #         end
    #     end
    #     for (j,(l2,m2,n2)) in enumerate(lmn_p)
    #         l2,m2,n2 = T.((l2,m2,n2))
    #         _coriolis_ts(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
    #     end
    # end

    is2, js2, aijs2 = rhs2(eachindex(lmn_t), lmn_p, lmn_t, Ω, ν)
    append!(is,is2)
    append!(js,js2)
    append!(aijs,aijs2)
    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

# function dist_ind(indices, np)
#     ni = indices÷np


function rhsd(N,m; Ω::T = 2.0, ν::T = 1.0) where T
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    # @show np
    nt = length(lmn_t)
    nu = np+nt 
    f1 = i -> rhs1(i,lmn_p, lmn_t, Ω, ν)
    dist_indices = getindex.(distribute(lmn_p).indices,1)
    ijaijs = pmap(f1, dist_indices)
    # @show ijaijs[1]
    is,js,aijs = vcat(getindex.(ijaijs,1)...),vcat(getindex.(ijaijs,2)...) , vcat(getindex.(ijaijs,3)...)

    f2 = i -> rhs2(i,lmn_p, lmn_t, Ω, ν)
    dist_indices = getindex.(distribute(lmn_t).indices, 1)
    ijaijs2 = pmap(f2, dist_indices)
    is2, js2, aijs2 = vcat(getindex.(ijaijs2,1)...),vcat(getindex.(ijaijs2,2)...) , vcat(getindex.(ijaijs2,3)...)

    
    append!(is,is2)
    append!(js,js2)
    append!(aijs,aijs2)
 
    RHS = sparse(is, js, aijs, nu, nu)
    return RHS

end

function rhst(N,m; Ω = 2.0, ν = 1.0)
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    @show np
    nt = length(lmn_t)
    nu = np+nt

    nth = Threads.nthreads()

    

    is, js, aijs = [Int[] for _=1:nth] , [Int[] for _=1:nth], [ComplexF64[] for _=1:nth]

    @inbounds @sync Threads.@threads for i in shuffle(eachindex(lmn_p))
        # (i,(l,m,n)) in enumerate(lmn_p)
        it = Threads.threadid()
        l,m,n = Float64.(lmn_p[i])
        push!(is[it],i)
        push!(js[it],i)
        push!(aijs[it],_coriolis_ss(l,m; Ω))
        for (j, (l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = Float64.((l2,m2,n2))
            if (l==l2) && (m==m2)
                _viscous_ss(is[it],js[it],aijs[it],i,j, l,m,n,n2; ν)
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            l2,m2,n2 = Float64.((l2,m2,n2))
            _coriolis_st(is[it],js[it],aijs[it], i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end

    @inbounds @sync Threads.@threads for i in shuffle(eachindex(lmn_t))
        it = Threads.threadid()
        l,m,n = Float64.(lmn_t[i])
        push!(is[it],i)
        push!(js[it],i)
        push!(aijs[it],_coriolis_tt(l,m; Ω))
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
                l2,m2,n2 = Float64.((l2,m2,n2))
                _viscous_tt(is[it],js[it],aijs[it], i+np,j+np, l,m,n,n2; ν)  
            end
        end
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = Float64.((l2,m2,n2))
            _coriolis_ts(is[it],js[it],aijs[it], i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end

    RHS = sparse(vcat(is...),vcat(js...),vcat(aijs...), nu, nu)
    return RHS

end



function rhs_coriolis(N,m; Ω = 2.0)
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], ComplexF64[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        l,m,n = Float64.((l,m,n))
        push!(is,i)
        push!(js,i)
        push!(aijs,_coriolis_ss(l,m; Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_t)
            l2,m2,n2 = Float64.((l2,m2,n2))
            _coriolis_st(is,js,aijs, i,j+np, l,l2,m,m2,n,n2; Ω)
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = Float64.((l,m,n))
        push!(is,i+np)
        push!(js,i+np)
        push!(aijs,_coriolis_tt(l,m; Ω))
        for (j,(l2,m2,n2)) in enumerate(lmn_p)
            l2,m2,n2 = Float64.((l2,m2,n2))
            _coriolis_ts(is,js,aijs, i+np,j, l,l2,m,m2,n,n2; Ω)
        end
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end

function rhs_viscosity(N,m; ν = 1.0)
    lmn_p = lmn_upol(N,m)
    lmn_t = lmn_utor(N,m)

    np = length(lmn_p)
    nt = length(lmn_t)
    nu = np+nt

    is,js,aijs = Int[], Int[], ComplexF64[]


    for (i,(l,m,n)) in enumerate(lmn_p)
        for (j, (l2,m2,n2)) in enumerate(lmn_p)
            if (l==l2) && (m==m2)
                l,m,n,n2 = Float64.((l,m,n,n2))
                _viscous_ss(is,js,aijs,i,j, l,m,n,n2; ν)
            end
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
                l,m,n,n2 = Float64.((l,m,n,n2))
                _viscous_tt(is,js,aijs, i+np,j+np, l,m,n,n2; ν)  
            end
        end
    end

    RHS = sparse(is,js,aijs, nu, nu)
    return RHS

end


end #module