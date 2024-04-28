module Poly

using SpecialFunctions
using Wigxjpf
using ForwardDiff
using DocStringExtensions

export wigner3j, wigner6j, wigner9j, adamgaunt, elsasser, jacobi, ylm, jacobis, ∂, p, _∂ll, D, innert, inners


#only for dev, initiate wigner symbols temporary working arrays.
function __wiginit(N)
    wig_table_init(2N, 9)
    wig_temp_init(2N)
end

function __wiginit_thread(N)
    wig_table_init(2N, 9)
    wig_thread_temp_init(2N)
end

#convenience for half-integer notation.
wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9) = wig9jj(2j1,2j2,2j3,2j4,2j5,2j6,2j7,2j8,2j9)
wigner6j(j1,j2,j3,j4,j5,j6) = wig6jj(2j1,2j2,2j3,2j4,2j5,2j6)
wigner3j(j1,j2,j3,j4,j5,j6) = wig3jj(2j1,2j2,2j3,2j4,2j5,2j6)

"""
$(TYPEDSIGNATURES)

Adam-Gaunt integral \$ A_{abc} = ...\$.
"""
@inline function adamgaunt(la,lb,lc,ma,mb,mc)::ComplexF64
    return (-1)^(mc)*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*wigner3j(Int(la), Int(lb), Int(lc), 0, 0, 0)*wigner3j(Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc))
end

@inline _Δ(la,lb,lc) = sqrt((la+lb+lc+2)*(la+lb+lc+4)/(4*(la+lb+lc+3)))*sqrt(complex((la+lb-lc+1)*(la-lb+lc+1)*(-la+lb+lc+1)))

"""
$(TYPEDSIGNATURES)

Elsasser variable \$ E_{abc} = ...\$.
"""
@inline function elsasser(la,lb,lc,ma,mb,mc)::ComplexF64
    return -(-1)^(mc)*im*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*_Δ(la,lb,lc)*wigner3j(Int(la)+1, Int(lb)+1, Int(lc)+1, 0, 0, 0)*wigner3j(Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc)) 
end

"""
$(TYPEDSIGNATURES)

Jacobi polynomial \$J_n^{(a,b)}(x)\$.
"""
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

"""
$(TYPEDSIGNATURES)

Jacobi polynomials \$J_n^{(a,b)}(2r^2-1)\$ and derivatives in \$r\$ up to third order for all degrees \$n \\in [0,N]\$, evaluated at a given grid `rgrid` of values in \$r\$.
"""
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

"""
$(TYPEDSIGNATURES)

Jacobi polynomials \$J_n^{(a,l+1/2)}(2r^2-1)\$ and derivatives in \$r\$ up to third order for all degrees \$n \\in [0,N]\$ and \$l \\in [1,N]\$. Default is `a=0.0`.
"""
function jacobis_l(N,r,a=0.0)
    return [jacobis(N,a,l+1/2,r) for l in 1:N] 
end

"""
$(TYPEDSIGNATURES)

Spherical harmonic in full norm, i.e. \$\\int Y_l^mY_i^j\\, \\sin(\\theta)\\,\\mathrm{d}\\theta\\mathrm{d}\\phi = \\delta_{li}\\delta_{mj}\$.
"""
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

"""
$(TYPEDSIGNATURES)

Spherical harmonic derivative in \$\\theta\$.
"""
function dylmdθ(l,m,θ,ϕ)
    return m*cot(θ)*ylm(l,m,θ,ϕ) + sqrt((l-m)*(l+m+1))*exp(-im*ϕ)*ylm(l,m+1,θ,ϕ)  
end
"""
$(TYPEDSIGNATURES)

Spherical harmonic derivative in \$\\phi\$.
"""
function dylmdϕ(l,m,θ,ϕ)
    return im*m*ylm(l,m,θ,ϕ)
end


"""
$(TYPEDSIGNATURES)

\$l(l+1)\$
"""
@inline p(l) = l*(l+1.0)


"""
$(TYPEDSIGNATURES)

Equation (25) in Ivers & Phillips (2008).
"""
function _∂ll(f,l,l1,r) 
    @assert l1 ∈ (l-1, l+1)
    if l1 == l-1
        return ∂(f,r) + (l+1)/r*f(r) 
    elseif l1 == l+1
        return ∂(f,r)-l/r*f(r)
    end
end

"""
$(TYPEDSIGNATURES)

Equation (xx) in Ivers & Phillips (2008).
"""
@inline D(f,l,r) = ∂(r->∂(f,r),r) + 2/r * ∂(f,r) - l*(l+1)/r^2 *f(r)

"""
$(TYPEDSIGNATURES)

Equation (xx) in Ivers & Phillips (2008).
"""
@inline function innert(t,t2, l, r) 
    return l*(l+1)*t(r)*t2(r)
end

"""
$(TYPEDSIGNATURES)

Equation (xx) in Ivers & Phillips (2008).
"""
@inline function inners(s,s2, l, r) 
    return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r*s(r),r)*∂(r->r*s2(r),r))/r^2
end

end #module