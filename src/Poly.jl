module Poly

using SpecialFunctions
using WignerSymbols
using ForwardDiff
using DocStringExtensions

export wigner3j, adamgaunt, elsasser, jacobi, ylm, ∂, p, _∂ll, D, dD, innert, inners, derivatives01, derivatives012, derivatives0123


"""
    adamgaunt(la,lb,lc,ma,mb,mc)

Adam-Gaunt integral [james_adams_1973](@citep), given by

```math
A_{abc} = \\oint\\int Y_aY_bY_c\\sin\\theta\\,\\mathrm{d}\\theta\\mathrm{d}\\phi,
```
where the spherical harmonics are abbreviated, so that ``Y_a = Y_{l_a}^{m_a}``.
"""
@inline function adamgaunt(la,lb,lc,ma,mb,mc)::ComplexF64
    return (-1)^(mc)*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*wigner3j(Float64,Int(la), Int(lb), Int(lc), 0, 0, 0)*wigner3j(Float64,Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc))
end

@inline _Δ(la,lb,lc) = sqrt((la+lb+lc+2)*(la+lb+lc+4)/(4*(la+lb+lc+3)))*sqrt(complex((la+lb-lc+1)*(la-lb+lc+1)*(-la+lb+lc+1)))

"""
    elsasser(la,lb,lc,ma,mb,mc)

Elsasser integral [james_adams_1973](@citep), given by

```math
E_{abc} = \\oint\\int Y_c\\left( \\frac{\\partial Y_a}{\\partial \\theta} \\frac{\\partial Y_b}{\\partial \\phi} - \\frac{\\partial Y_a}{\\partial \\phi}\\frac{\\partial Y_b}{\\partial \\theta} \\right)\\,\\mathrm{d}\\theta\\mathrm{d}\\phi,
```
where the spherical harmonics are abbreviated, so that ``Y_a = Y_{l_a}^{m_a}``.
"""
@inline function elsasser(la,lb,lc,ma,mb,mc)::ComplexF64
    return -(-1)^(mc)*im*sqrt((2la + 1)*(2lb + 1)*(2lc + 1)/4π)*_Δ(la,lb,lc)*wigner3j(Float64,Int(la)+1, Int(lb)+1, Int(lc)+1, 0, 0, 0)*wigner3j(Float64,Int(la),Int(lb),Int(lc),Int(ma),Int(mb),-Int(mc)) 
end

"""
    jacobi(n,a,b,x)

Jacobi polynomial 

```math
J_n^{(a,b)}(x)
```
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

using ForwardDiff: Dual, Tag, value, partials

function derivatives01(f::F, x::T) where {F, T<:Real}
    # Outer perturbation (carries f')
    Touter = typeof(Tag(f, T))
    x1 = Dual{Touter}(x, one(x))            # x + δ

    y = f(x1)                               # <-- the only call to f

    f0 = value(y)
    f1 = partials(y, 1)
    return f0, f1
end

function derivatives012(f::F, x::T) where {F, T<:Real}
    # Outer perturbation (carries f')
    Touter = typeof(Tag(f, T))
    x1 = Dual{Touter}(x, one(x))            # x + δ

    # Inner perturbation, nested over the outer dual (carries f'')
    Tinner = typeof(Tag(f, typeof(x1)))
    x2 = Dual{Tinner}(x1, one(x1))          # (x + δ) + ε

    y = f(x2)                               # <-- the only call to f

    f0 = value(value(y))
    f1 = value(partials(y, 1))
    f2 = partials(partials(y, 1), 1)
    return f0, f1, f2
end

function derivatives0123(f::F, x::T) where {F, T<:Real}
    T1 = typeof(Tag(f, T))
    x1 = Dual{T1}(x, one(x))            # x + ε₁

    T2 = typeof(Tag(f, typeof(x1)))
    x2 = Dual{T2}(x1, one(x1))          # x + ε₁ + ε₂

    T3 = typeof(Tag(f, typeof(x2)))
    x3 = Dual{T3}(x2, one(x2))          # x + ε₁ + ε₂ + ε₃

    y = f(x3)

    f0 = value(value(value(y)))
    f1 = value(value(partials(y, 1)))
    f2 = value(partials(partials(y, 1), 1))
    f3 = partials(partials(partials(y, 1), 1), 1)
    return f0, f1, f2, f3
end

"""
$(TYPEDSIGNATURES)

Spherical harmonic ``Y_l^m`` in full norm, i.e. 
```math
\\int Y_l^mY_i^j\\, \\sin(\\theta)\\,\\mathrm{d}\\theta\\mathrm{d}\\phi = \\delta_{li}\\delta_{mj}
```
where ``\\theta`` is the colatitude and ``\\phi`` the azimuthal angle.
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

Derivative of spherical harmonic ``Y_l^m`` in ``\\theta``.
"""
function dylmdθ(l,m,θ,ϕ)
    return m*cot(θ)*ylm(l,m,θ,ϕ) + sqrt((l-m)*(l+m+1))*exp(-im*ϕ)*ylm(l,m+1,θ,ϕ)  
end
"""
$(TYPEDSIGNATURES)

Derivative of spherical harmonic ``Y_l^m`` in ``\\phi``.
"""
function dylmdϕ(l,m,θ,ϕ)
    return im*m*ylm(l,m,θ,ϕ)
end


"""
$(TYPEDSIGNATURES)

``p(l) = l(l+1)`` following notation of [ivers_scalar_2008](@citet)
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
    _f, _df = derivatives01(f,r)
    if l1 == l-1
        return _df + (l+1)/r*_f
    elseif l1 == l+1
        return _df-l/r*_f
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
@inline function D(f,l,r)
    _f, _df, _d2f = derivatives012(f,r)
    return D(_f, _df, _d2f,l,r)
end
# @inline D(f,l,r) = ∂(r->∂(f,r),r) + 2/r * ∂(f,r) - l*(l+1)/r^2 *f(r)
@inline function D(_f, _df, _d2f,l,r)
    return _d2f + 2/r*_df - l*(l+1)/r^2*_f
end

@inline function dD(_f, _df, _d2f, _d3f, l, r)
    return _d3f -2/r^2*_df + 2/r*_d2f + 2l*(l+1)/r^3*_f - l*(l+1)/r^2*_df
end

"""
$(TYPEDSIGNATURES)

```math
r\\rightarrow l(l+1) t(r) t_2(r)
```

Radial function to be integrated in radius when computing the inner product of two toroidal vectors.
"""
@inline function innert(t::T1,t2::T2, l::Int, r::Tr) where {T1,T2,Tr}
    return l*(l+1)*t(r)*t2(r)
end

"""
$(TYPEDSIGNATURES)

```math
r\\rightarrow \\frac{l(l+1)}{r^2}\\left( l(l+1)s(r) s_2(r) + \\frac{\\partial r s(r)}{\\partial r}\\frac{\\partial r s_2(r)}{\\partial r}\\right)
```

Radial function to be integrated in radius when computing the inner product of two poloidal vectors.
"""
@inline function inners(s,s2, l, r) 
    _s, _ds = derivatives01(s,r)
    _s2, _ds2 = derivatives01(s2,r)
    return l*(l+1)*(_s*_s2*l*(l+1)+(r*_ds+_s)*(r*_ds2+_s2))/r^2
end

end #module
