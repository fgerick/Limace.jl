@inline Base.@pure function jacobi(n,a,b,x)
    ox = one(x)
    zx = zero(x)
    if n==0
        return ox
    elseif n==1
        return ox/2 * (a - b + (a + b + 2)*x)
    end

    p0 = ox
    p1 = ox/2 * (a - b + (a + b + 2)*x)
    p2 = zx;

    for i = 1:(n-1)
        a1 = 2*(i+1)*(i+a+b+1)*(2*i+a+b);
        a2 = (2*i+a+b+1)*(a*a-b*b);
        a3 = (2*i+a+b)*(2*i+a+b+1)*(2*i+a+b+2);
        a4 = 2*(i+a)*(i+b)*(2*i+a+b+2);
        p2 = ox/a1*( (a2 + a3*x)*p1 - a4*p0);

        p0 = p1
        p1 = p2
    end

    return p2
end

@inline djacobidr(n,a,b,r) = 2(1+a+b+n)*r*jacobi(n-1,1+a,1+b,2r^2-1)
@inline djacobid2r(n,a,b,r) = 2(1+a+b+n)*jacobi(n-1,1+a,1+b,2r^2-1) + 2*(1+a+b+n)*r*djacobidr(n-1,a+1,b+1,r) 
@inline djacobid3r(n,a,b,r) = 2(1+a+b+n)*djacobidr(n-1,1+a,1+b,r) + 2*(1+a+b+n)*djacobidr(n-1,a+1,b+1,r) + 2*(1+a+b+n)*r*djacobid2r(n-1,a+1,b+1,r)

@inline function djacobi(n,a,b,x,k)
    exp(loggamma(a+b+n+1+k)-loggamma(a+b+n+1))/2^k*jacobi(n-k,a+k,b+k,x)
end

function jacobis(N,a,b,r,kmax)
    nr = length(r)
    js = zeros(eltype(r),nr,N,kmax+1)
    # for n in 1:N
    #     js[:,n,1] .= jacobi.(n,a,b,r)
    # end
    for k in 0:kmax, n in 1:N
        js[:,n, k+1] .= djacobi.(n,a,b,r,k)
    end
    return js
end

