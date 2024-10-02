module Poly

using SpecialFunctions
using WignerSymbols
using ForwardDiff
using DocStringExtensions

export wigner3j, adamgaunt, elsasser, jacobi, ylm, jacobis, ∂, p, _∂ll, D, innert, inners, _∂ll_m1, _∂ll_p1


"""
$(TYPEDSIGNATURES)

Adam-Gaunt integral \$ A_{abc} = \\oint\\int Y_iY_jY_k\\sin\\theta\\,\\mathrm{d}\\theta\\mathrm{d}\\phi\$.
"""
@inline function adamgaunt(la,lb,lc,ma,mb,mc)::ComplexF64
    return (-1)^(mc)*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*wigner3j(Float64,Int(la), Int(lb), Int(lc), 0, 0, 0)*wigner3j(Float64,Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc))
end

@inline _Δ(la,lb,lc) = sqrt((la+lb+lc+2)*(la+lb+lc+4)/(4*(la+lb+lc+3)))*sqrt(complex((la+lb-lc+1)*(la-lb+lc+1)*(-la+lb+lc+1)))

"""
$(TYPEDSIGNATURES)

Elsasser variable \$ E_{abc} = \\oint\\int Y_k\\left( \\frac{\\partial Y_i}{\\partial \\theta} \\frac{\\partial Y_j}{\\partial \\phi} - \\frac{\\partial Y_i}{\\partial \\phi}\\frac{\\partial Y_j}{\\partial \\theta} \\right)\\,\\mathrm{d}\\theta\\mathrm{d}\\phi\$\$.
"""
@inline function elsasser(la,lb,lc,ma,mb,mc)::ComplexF64
    return -(-1)^(mc)*im*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*_Δ(la,lb,lc)*wigner3j(Float64,Int(la)+1, Int(lb)+1, Int(lc)+1, 0, 0, 0)*wigner3j(Float64,Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc)) 
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

Derivative of spherical harmonic in \$\\theta\$.
"""
function dylmdθ(l,m,θ,ϕ)
    return m*cot(θ)*ylm(l,m,θ,ϕ) + sqrt((l-m)*(l+m+1))*exp(-im*ϕ)*ylm(l,m+1,θ,ϕ)  
end
"""
$(TYPEDSIGNATURES)

Derivative of spherical harmonic in \$\\phi\$.
"""
function dylmdϕ(l,m,θ,ϕ)
    return im*m*ylm(l,m,θ,ϕ)
end


"""
$(TYPEDSIGNATURES)

\$p(l) = l(l+1)\$ following [ivers_scalar_2008](@citet)
"""
@inline p(l) = l*(l+1.0)


"""
$(TYPEDSIGNATURES)

```math
\\partial_l^{l_1} = \\begin{cases}
\\frac{\\partial f}{\\partial r} + \\frac{l+1}{r}f \\quad \\mathrm{if}\\, l_1 = l-1\\\\
\\frac{\\partial f}{\\partial r} - \\frac{l}{r}f \\quad \\mathrm{if}\\, l_1 = l+1
\\end{cases}
```

Equation (25) in [ivers_scalar_2008](@citet).
"""
@inline function _∂ll(f,l,l1,r)
    # @assert l1 ∈ (l-1, l+1)
    if l1 == l-1
        return ∂(f,r) + (l+1)/r*f(r) 
    elseif l1 == l+1
        return ∂(f,r)-l/r*f(r)
    end
    return 0.0
end

"""
$(TYPEDSIGNATURES)

```math
D_l(f) = \\frac{\\partial^2 f}{\\partial r^2} +\\frac{2}{r}\\frac{\\partial f}{\\partial r} - \\frac{l(l+1)}{r^2}f
```

Below equation (31) in [ivers_scalar_2008](@citet).
"""
@inline D(f,l,r) = ∂(r->∂(f,r),r) + 2/r * ∂(f,r) - l*(l+1)/r^2 *f(r)

"""
$(TYPEDSIGNATURES)

```math
l(l+1)t t_2
```
"""
@inline function innert(t::T1,t2::T2, l::Int, r::Tr) where {T1,T2,Tr}
    return l*(l+1)*t(r)*t2(r)
end

"""
$(TYPEDSIGNATURES)

```math
\\frac{l(l+1)}{r^2}\\left( l(l+1)s s_2 + \\frac{\\partial r s}{\\partial r}\\frac{\\partial r s_2}{\\partial r}\\right)
```
"""
@inline function inners(s,s2, l, r) 
    return l*(l+1)*(s(r)*s2(r)*l*(l+1)+∂(r->r*s(r),r)*∂(r->r*s2(r),r))/r^2
end

end #module