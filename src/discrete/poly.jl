wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) = wig9jj(2j1,2j2,2j3,2j4,2j5,2j6,2j7,2j8,2j9)
wigner6j(j1,j2,j3,j4,j5,j6) = wig6jj(2j1,2j2,2j3,2j4,2j5,2j6)
wigner3j(j1,j2,j3,j4,j5,j6) = wig3jj(2j1,2j2,2j3,2j4,2j5,2j6)


function adamgaunt(la,lb,lc,ma,mb,mc)
    return -sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*wigner3j(la, lb, lc, 0, 0, 0)*wigner3j(la,lb,lc,ma,mb,-mc)
end

Δ(la,lb,lc) = sqrt((la+lb+lc+2)*(la+lb+lc+4)/(4*(la+lb+lc+3)))*sqrt((la+lb-lc+1)*(la-lb+lc+1)*(-la+lb+lc+1))

function elsasser(la,lb,lc,ma,mb,mc)
    return im*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*Δ(la,lb,lc)*wigner3j(la+1, lb+1, lc+1, 0, 0, 0)*wigner3j(la,lb,lc,ma,mb,-mc) 
end

@inline function jacobi(n,a,b,x)
    ox = one(x)
    zx = zero(x)
    if n==0
        return ox
    elseif n==1
        return ox/2 * (a - b + (a + b + 2)*x)
    elseif n<0 #convenience
        return zx
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


const ∂ =  ForwardDiff.derivative