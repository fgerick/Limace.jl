@inline function t_chen(l,m,n,r) 
    fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n)))
    return fac * r^l * (jacobi(n,0,l+1/2, 2r^2-1) - jacobi(n-1,0,l+1/2,2r^2-1)) 
end

#2.38,2.39 poloidal scalar
@inline function s_chen(l,m,n,r)
    c1 = 2l+4n+1
    c2 = -2(2l+4n+3)
    c3 = 2l+4n+5
    fac = 1/(sqrt(2l*(1 + l)*(1 + 2*l + 4*n)*(3 + 2*l + 4*n)*(5 + 2*l + 4*n)))
    return fac*r^l*(c1*jacobi(n+1,0,l+1/2,2r^2-1) + c2*jacobi(n,0,l+1/2,2r^2-1) + c3*jacobi(n-1,0,l+1/2,2r^2-1)  ) 
end

@inline function s_mf(l, m, n, r)
    fac = (n==1) ? 1/sqrt(l*(1 + l)*(5 + 2*l)*(6 + l*(11 + 6*l))) : 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    return fac * r^l * ( (2*l + 4*n - 3) * jacobi(n,0,l+1/2,2*r^2-1) - 2*(2*l + 4*n - 1)*jacobi(n-1,0,l+1/2,2*r^2-1) +(2*l + 4*n + 1)*jacobi(n-2,0,l+1/2,2*r^2-1))
end

const t_mf = t_chen

# @inline function t_mf(l, m, n, r) 
#     fac = 1/sqrt(l*(1 + l)*(1/(-1 + 2*l + 4*n) + 1/(3 + 2*l + 4*n)))
#     return fac * r^l * (jacobi(n,0,l+1/2,2*r^2-1) - jacobi(n-1,0,l+1/2,2*r^2-1))
# end


# @inline function t_visc(l,m,n,r) 
#     fac = 1/2 * sqrt((3+2l+2n)*(5+2l+2n)*(7+2l+4n))/sqrt(l*(l+1)*(n+1)*(n+2))
#     return (1-r^2)*r^l*jacobi(n,2,l+1/2, 2r^2-1)*fac
# end

# @inline function s_visc(l,m,n,r) 
#     fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
#     return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
# end

# @inline function drs_visc(l,m,n,r) 
#     s = s_visc(l,m,n,r)
#     fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
#     return -2r^2*s/(1-r^2) + (l+1)*s + (1-r^2)*r^(l+1)*djacobidr(n,1,l+1/2,r)*fac
# end

# @inline function ds_visc(l,m,n,r) 
#     s = s_visc(l,m,n,r)
#     fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
#     return -2r*s/(1-r^2) + l*s/r + (1-r^2)*r^l*djacobidr(n,1,l+1/2,r)*fac
# end

# @inline function dt_visc(l,m,n,r) 
#     t = t_visc(l,m,n,r)
#     fac = 1/2 * sqrt((3+2l+2n)*(5+2l+2n)*(7+2l+4n))/sqrt(l*(l+1)*(n+1)*(n+2))
#     return -2r*t/(1-r^2) + l*t/r + (1-r^2)*r^l*djacobidr(n,2,l+1/2,r)*fac
# end