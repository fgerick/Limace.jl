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
    fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    return fac * r^l * ( (2*l + 4*n - 3) * jacobi(n,0,l+1/2,2*r^2-1) - 2*(2*l + 4*n - 1)*jacobi(n-1,0,l+1/2,2*r^2-1) +(2*l + 4*n + 1)*jacobi(n-2,0,l+1/2,2*r^2-1))
end

@inline function ds_mf_dr(l,m,n,r)
    fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    t1 = s_mf(l,m,n,r)/r*l 
    t2 = fac * r^l * ( (2*l + 4*n - 3) * djacobidr(n,0,l+1/2,r) - 2*(2*l + 4*n - 1)*djacobidr(n-1,0,l+1/2,r) +(2*l + 4*n + 1)*djacobidr(n-2,0,l+1/2,r))
    return t1 + t2
end

@inline function ds_mf_d2r(l,m,n,r)
    fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))
    t1 = -s_mf(l,m,n,r)/r^2*l
    t1 += ds_mf_dr(l,m,n,r)/r*l
    t2 = fac * r^l * ( (2*l + 4*n - 3) * djacobid2r(n,0,l+1/2,r) - 2*(2*l + 4*n - 1)*djacobid2r(n-1,0,l+1/2,r) +(2*l + 4*n + 1)*djacobid2r(n-2,0,l+1/2,r))
    t2 += fac * r^(l-1) * l * ( (2*l + 4*n - 3) * djacobidr(n,0,l+1/2,r) - 2*(2*l + 4*n - 1)*djacobidr(n-1,0,l+1/2,r) +(2*l + 4*n + 1)*djacobidr(n-2,0,l+1/2,r)) 
    return t1 + t2
end

# @inline function ds_mf_d3r(l,m,n,r)
#     fac = 1/(sqrt(2l*(1 + l)*(-3 + 2*l + 4*n)*(-1 + 2*l + 4*n)*(1 + 2*l + 4*n)))

#     t1 = s_mf(l,m,n,r)/r^3*l*2
#     t1 += -ds_mf_dr(l,m,n,r)/r^2*l

#     t1 -=  ds_mf_dr(l,m,n,r)/r^2*l
#     t1 +=  ds_mf_d2r(l,m,n,r)/r*l
     
#     t2 = fac * r^(l-2) * l*(l-1) * ( (2*l + 4*n - 3) * djacobid3r(n,0,l+1/2,r) - 2*(2*l + 4*n - 1)*djacobid3r(n-1,0,l+1/2,r) +(2*l + 4*n + 1)*djacobid3r(n-2,0,l+1/2,r))
#     return t1 + t2
# end

const t_mf = t_chen


@inline function s_mf1(l, m, n, r)
    c0 = -2n^2*(l + 1) - n*(l + 1)*(2l - 1) - l*(2l + 1)
    c1 = (2*(l + 1)*n^2 + (2l + 3)*(l + 1)*n + (2l + 1)^2)
    c2 = (4n*l + l*(2l+1))
    fac = 1/sqrt(l*(l+1.0)*(1.0+2l+4n)*((1.0+2l)^2 + (l+1)*(3+2l)*n + 2(1.0+l)*n^2) * (2l^2*(n+1.0) + n*(2n-1.0) + l*(2n^2+n+1.0)) )
    return (fac*c0 * jacobi(n, 0.0, 1 / 2 + l, -1 + 2 * r ^ 2) + fac*c1 * jacobi(-1 + n, 0, 1 / 2 + l, -1 + 2 * r ^ 2) + fac*c2)*r^l
end


@inline function t_mf1(l, m, n, r) 
    fac = 1/2*sqrt(9.0+8l^3+36n+44n^2+16n^3+4l^2*(7.0+8n)+l*(30.0+72n+40.0n^2))/sqrt(l*(l+1.0)*n*(n+1))
    return (1 - r^2) * jacobi(n - 1, 2.0, l + 1/2, 2r^2 - 1) * r^l * fac
end

#inviscid velocities
#https://homepages.see.leeds.ac.uk/~earpwl/Galerkin/Galerkin.html (5.1, 5.6)
function t_in(l,m,n,r) 
    fac = sqrt(3+2l+4n)/sqrt(l*(l+1))
    return r^l*jacobi(n,0,l+1/2, 2r^2-1)*fac
end


function s_in(l,m,n,r) 
    fac = sqrt(5+2l+4n)/sqrt(4l*(l+1)*(n+1)^2)
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end


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