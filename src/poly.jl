module Poly

using SpecialFunctions
using Wigxjpf
using ForwardDiff

export wigner3j, wigner6j, wigner9j, adamgaunt, elsasser, jacobi, ylm, jacobis, ∂

wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) = wig9jj(2j1,2j2,2j3,2j4,2j5,2j6,2j7,2j8,2j9)
wigner6j(j1,j2,j3,j4,j5,j6) = wig6jj(2j1,2j2,2j3,2j4,2j5,2j6)
wigner3j(j1,j2,j3,j4,j5,j6) = wig3jj(2j1,2j2,2j3,2j4,2j5,2j6)


@inline function adamgaunt(la,lb,lc,ma,mb,mc)
    return (-1)^(mc)*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*wigner3j(Int(la), Int(lb), Int(lc), 0, 0, 0)*wigner3j(Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc))
end

@inline Δ(la,lb,lc) = sqrt((la+lb+lc+2)*(la+lb+lc+4)/(4*(la+lb+lc+3)))*sqrt(complex((la+lb-lc+1)*(la-lb+lc+1)*(-la+lb+lc+1)))


@inline function elsasser(la,lb,lc,ma,mb,mc)
    return -(-1)^(mc)*im*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*Δ(la,lb,lc)*wigner3j(Int(la)+1, Int(lb)+1, Int(lc)+1, 0, 0, 0)*wigner3j(Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc)) 
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
        _2iab = 2i+a+b
        a1 = 2*(i+1)*(i+a+b+1)*_2iab
        a2 = (_2iab+1)*(a*a-b*b)
        a3 = _2iab*(_2iab+1)*(_2iab+2)
        a4 = 2*(i+a)*(i+b)*(_2iab+2)
        p2 = ox/a1*( (a2 + a3*x)*p1 - a4*p0)

        p0 = p1
        p1 = p2
    end

    return p2
end


const ∂ =  ForwardDiff.derivative

#derivatives of Jacobi polynomials ∂/∂r(jacobi(n,a,b, x = 2r^2-1)).
@inline djacobidr(n,a,b,r) = 2(1+a+b+n)*r*jacobi(n-1,1+a,1+b,2r^2-1)
@inline djacobid2r(n,a,b,r) = 2(1+a+b+n)*jacobi(n-1,1+a,1+b,2r^2-1) + 2*(1+a+b+n)*r*djacobidr(n-1,a+1,b+1,r) 
@inline djacobid3r(n,a,b,r) = 2(1+a+b+n)*djacobidr(n-1,1+a,1+b,r) + 2*(1+a+b+n)*djacobidr(n-1,a+1,b+1,r) + 2*(1+a+b+n)*r*djacobid2r(n-1,a+1,b+1,r)

@inline function djacobi(n,a,b,x,k)
    exp(loggamma(a+b+n+1+k)-loggamma(a+b+n+1))/2^k*jacobi(n-k,a+k,b+k,x)
end

function jacobis(N,a,b,rgrid)
    nr = length(rgrid)
    js = zeros(eltype(rgrid),4,N+1,nr) 
    @inbounds for n in 0:N, i=1:nr
        r = rgrid[i]
        js[1,n+1,i] = jacobi(n,a,b, 2r^2-1)
        js[2,n+1,i] = djacobidr(n,a,b,r)
        js[3,n+1,i] = djacobid2r(n,a,b,r)
        js[4,n+1,i] = djacobid3r(n,a,b,r)
    end
    return js   
end

function jacobis_l(N,r,a=0.0)
    return [jacobis(N,a,l+1/2,r) for l in 1:N] 
end

# radial derivative ∂/∂r(r^l J) with J a jacobi polynomial
@inline function d_rlJ(l, r, rl, J, dJ)
    # rl = r^l
	return rl*(l/r*J + dJ)
end

@inline function d2_rlJ(l, r, rl, J, dJ, d2J)
    # rl = r^l
	return rl*(l*(l-1)/r^2*J + 2l/r*dJ + d2J)
end

@inline function d3_rlJ(l, r, rl, J, dJ, d2J, d3J)
    # rl = r^l
	return rl*(l*(l-1)*(l-2)/r^3*J + 3l*(l-1)/r^2*dJ + 3l/r*d2J + d3J)
end

function ylm(ℓ::Int, m::Int, θ, φ) #norm -> ∫YₗᵐYᵢʲsin(θ)dθdϕ = δₗᵢδₘⱼ
    if ℓ<abs(m)
        return zero(complex(typeof(θ)))
    else
        m̃ = abs(m)
        a =  exp((loggamma(ℓ+m̃+1)+loggamma(ℓ-m̃+1)-2loggamma(ℓ+1))/2) *sqrt(2ℓ+1)/sqrt(4π)
        if m<0
            a*=(-1)^m
        end
        return a * exp(im*m*φ) * (-sin(θ/2) * cos(θ/2))^m̃ * jacobi(ℓ-m̃,m̃,m̃,cos(θ))
    end
end

function dylmdθ(l,m,θ,ϕ)
    return m*cot(θ)*ylm(l,m,θ,ϕ) + sqrt((l-m)*(l+m+1))*exp(-im*ϕ)*ylm(l,m+1,θ,ϕ)  
end

function dylmdϕ(l,m,θ,ϕ)
    return im*m*ylm(l,m,θ,ϕ)
end


function poloidal_discretize(s,l,m,n,r,θ,ϕ)
    ur = l*(l+1)*s(l,m,n,r)*ylm(l,m,θ,ϕ)/r
    uθ = 1/r*∂(r->s(l,m,n,r)*r,r)*dylmdθ(l,m,θ,ϕ)
    uϕ = 1/(r*sin(θ))*∂(r->s(l,m,n,r)*r,r)*dylmdϕ(l,m,θ,ϕ)
    return (ur,uθ,uϕ)
end

function toroidal_discretize(t,l,m,n,r,θ,ϕ)
    ur = 0.0
    uθ = 1/sin(θ)*t(l,m,n,r)*dylmdϕ(l,m,θ,ϕ)
    uϕ = -t(l,m,n,r)*dylmdθ(l,m,θ,ϕ)

    return (ur,uθ,uϕ)
end

end #module