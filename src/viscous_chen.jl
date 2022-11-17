module ChenBasis

using SparseArrays

lmn_upol(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:(N-1) for n in 1:((N-l+1)÷2+1) if abs(m)<=l] 

lmn_utor(N, ms = 0:N) = [(l,m,n) for m in ms for l in 1:N for n in 1:((N-l)÷2+1) if abs(m)<=l] 

function _inner_tt(is,js,aijs, i,j, l,n,n2; ν = 1.0) 
    
   
    if n==n2 
        aij = 1.0
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

# 1/(2.*sqrt(1 + 3/(2 + 4*l + 8*n) - 3/(18 + 4*l + 8*n)))

function _inner_ss(is,js,aijs, i,j, l,n,n2; ν = 1.0) 
    
   
    if n==n2 
        aij = 1.0
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

function lhs(N,m; Ω::T = 1.0) where T
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
                _inner_ss(is,js,aijs,i,j, l,n,n2)
            end
        end
    end

    for (i,(l,m,n)) in enumerate(lmn_t)
        l,m,n = T.((l,m,n))
        for (j, (l2,m2,n2)) in enumerate(lmn_t)
            if (l==l2) && (m==m2)
                l2,m2,n2 = T.((l2,m2,n2))
                _inner_tt(is,js,aijs, i+np,j+np, l,n,n2)  
            end
        end
    end

    LHS = sparse(is,js,aijs, nu, nu)
    return LHS

end


function _coriolis_tt(is,js,aijs, i,j, l,m,n,n2; Ω = 2) 
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
function _coriolis_ss(is,js,aijs, i,j, l,m,n,n2; Ω = 2) 
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


function _coriolis_ts(is,js,aijs, i,j, l,l2,m,m2,n,n2; Ω = 2) 
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

function _coriolis_st(is,js,aijs, i,j, l2,l,m2,m,n2,n; Ω = 2) 
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


function _viscous_tt(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
    if n==n2 
        aij = -((-1 + 2*l + 4*n)*(3 + 2*l + 4*n))/2
        push!(is,i)
        push!(js,j)
        push!(aijs,ν*aij)
    end
        
    return nothing
end


function _viscous_ss(is,js,aijs, i,j, l,m,n,n2; ν = 1.0) 
    
   
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

end