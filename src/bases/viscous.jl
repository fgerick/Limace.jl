@inline function t_visc(l,m,n,r) 
    fac = 1/2 * sqrt((3+2l+2n)*(5+2l+2n)*(7+2l+4n))/sqrt(l*(l+1)*(n+1)*(n+2))
    return (1-r^2)*r^l*jacobi(n,2,l+1/2, 2r^2-1)*fac
end

@inline function s_visc(l,m,n,r) 
    fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
    return (1-r^2)*r^l*jacobi(n,1,l+1/2, 2r^2-1)*fac
end

@inline function drs_visc(l,m,n,r) 
    s = s_visc(l,m,n,r)
    fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
    return -2r^2*s/(1-r^2) + (l+1)*s + (1-r^2)*r^(l+1)*djacobidr(n,1,l+1/2,r)*fac
end

@inline function ds_visc(l,m,n,r) 
    s = s_visc(l,m,n,r)
    fac = sqrt(5+2l+4n)/sqrt(l*(l+1))/(2*(n+1))
    return -2r*s/(1-r^2) + l*s/r + (1-r^2)*r^l*djacobidr(n,1,l+1/2,r)*fac
end

@inline function dt_visc(l,m,n,r) 
    t = t_visc(l,m,n,r)
    fac = 1/2 * sqrt((3+2l+2n)*(5+2l+2n)*(7+2l+4n))/sqrt(l*(l+1)*(n+1)*(n+2))
    return -2r*t/(1-r^2) + l*t/r + (1-r^2)*r^l*djacobidr(n,2,l+1/2,r)*fac
end