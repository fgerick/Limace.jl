
#inviscid poloidal scalar

Base.@propagate_inbounds s_in_fac(l, n) = sqrt(5 + 2l + 4n) / sqrt(4l * (l + 1) * (n + 1)^2)

Base.@propagate_inbounds function s_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
    rl = rls[l][i]
    fac = s_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    return (1 - r^2) * rl * J * fac
end

Base.@propagate_inbounds function d_s_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = s_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    rlJ = rl * J
    return (-2r * rlJ + (1 - r^2) * d_rlJ(l, r, rl, J, dJ)) * fac
end

Base.@propagate_inbounds function d2_s_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = s_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    rlJ = rl * J
    return (-2rlJ - 4r * d_rlJ(l, r, rl, J, dJ) + (1 - r^2) * d2_rlJ(l, r, rl, J, dJ, d2J)) * fac
end

Base.@propagate_inbounds function d3_s_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = s_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    d3J = js[l][4, n+1, i]
    return (-6d_rlJ(l, r, rl, J, dJ) - 6r * d2_rlJ(l, r, rl, J, dJ, d2J) + (1 - r^2) * d3_rlJ(l, r, rl, J, dJ, d2J, d3J)) * fac
end


#inviscid toroidal scalar

Base.@propagate_inbounds t_in_fac(l, n) = sqrt(3 + 2l + 4n) / sqrt(l * (l + 1))
Base.@propagate_inbounds function t_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    # r = rgrid[i]
	rl = rls[l][i]
    J = js[l][1, n+1, i]
    fac = t_in_fac(T(l), T(n))
    return rl * J * fac
end


Base.@propagate_inbounds function d_t_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    return d_rlJ(l, r, rl, J, dJ) * fac
end

Base.@propagate_inbounds function d2_t_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    return d2_rlJ(l, r, rl, J, dJ, d2J) * fac
end

Base.@propagate_inbounds function d3_t_in_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_in_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    d3J = js[l][4, n+1, i]
    return d3_rlJ(l, r, rl, J, dJ, d2J, d3J) * fac
end

#viscous toroidal scalar

Base.@propagate_inbounds t_chen_fac(l, n) = 1 / sqrt(l * (1 + l) * (1 / (-1 + 2 * l + 4 * n) + 1 / (3 + 2 * l + 4 * n)))

Base.@propagate_inbounds function t_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    fac = t_chen_fac(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    Jn = js[l][1, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    return fac * rl * (Jn - Jnm1)
end


Base.@propagate_inbounds function d_t_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_chen_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]

    return (d_rlJ(l, r, rl, J, dJ) - d_rlJ(l, r, rl, Jnm1, dJnm1)) * fac
end

Base.@propagate_inbounds function d2_t_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_chen_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    return (d2_rlJ(l, r, rl, J, dJ, d2J) - d2_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1)) * fac
end

Base.@propagate_inbounds function d3_t_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    fac = t_chen_fac(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    d3J = js[l][4, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    d3Jnm1 = (n - 1 < 0) ? zero(T) : js[l][4, n, i]

    return (d3_rlJ(l, r, rl, J, dJ, d2J, d3J) - d3_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1, d3Jnm1)) * fac
end

#viscous poloidal scalar
Base.@propagate_inbounds function s_chen_fac_cs(l, n)
    c1 = 2l + 4n + 1
    c2 = -2(2l + 4n + 3)
    c3 = 2l + 4n + 5
    fac = 1 / (sqrt(2l * (1 + l) * (1 + 2 * l + 4 * n) * (3 + 2 * l + 4 * n) * (5 + 2 * l + 4 * n)))
    return c1, c2, c3, fac
end

Base.@propagate_inbounds function s_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    c1, c2, c3, fac = s_chen_fac_cs(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    Jnp1 = js[l][1, n+2, i]
    J = js[l][1, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]

    return fac * rl * (c1 * Jnp1 + c2 * J + c3 * Jnm1)
end

Base.@propagate_inbounds function d_s_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    c1, c2, c3, fac = s_chen_fac_cs(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    Jnp1 = (n - 1 < 0) ? zero(T) : js[l][1, n+2, i]
    dJnp1 = (n - 1 < 0) ? zero(T) : js[l][2, n+2, i]

    return (c1 * d_rlJ(l, r, rl, Jnp1, dJnp1) + c2 * d_rlJ(l, r, rl, J, dJ) + c3 * d_rlJ(l, r, rl, Jnm1, dJnm1)) * fac
end

Base.@propagate_inbounds function d2_s_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    c1, c2, c3, fac = s_chen_fac_cs(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    Jnp1 = (n - 1 < 0) ? zero(T) : js[l][1, n+2, i]
    dJnp1 = (n - 1 < 0) ? zero(T) : js[l][2, n+2, i]
    d2Jnp1 = (n - 1 < 0) ? zero(T) : js[l][3, n+2, i]

    return (c1 * d2_rlJ(l, r, rl, Jnp1, dJnp1, d2Jnp1) + c2 * d2_rlJ(l, r, rl, J, dJ, d2J) + c3 * d2_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1)) * fac
end

Base.@propagate_inbounds function d3_s_chen_pre(js, rls, l, m, n, rgrid::Vector{T}, i) where {T}
    r = rgrid[i]
	rl = rls[l][i]
    c1, c2, c3, fac = s_chen_fac_cs(T(l), T(n))
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    d3J = js[l][4, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    d3Jnm1 = (n - 1 < 0) ? zero(T) : js[l][4, n, i]
    Jnp1 = (n - 1 < 0) ? zero(T) : js[l][1, n+2, i]
    dJnp1 = (n - 1 < 0) ? zero(T) : js[l][2, n+2, i]
    d2Jnp1 = (n - 1 < 0) ? zero(T) : js[l][3, n+2, i]
    d3Jnp1 = (n - 1 < 0) ? zero(T) : js[l][4, n+2, i]

    return (c1 * d3_rlJ(l, r, rl, Jnp1, dJnp1, d2Jnp1, d3Jnp1) + c2 * d3_rlJ(l, r, rl, J, dJ, d2J, d3J) + c3 * d3_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1, d3Jnm1)) * fac
end


#insulating mf toroidal scalar
const t_mf_pre = t_chen_pre
const d_t_mf_pre = d_t_chen_pre
const d2_t_mf_pre = d2_t_chen_pre
const d3_t_mf_pre = d3_t_chen_pre

#insulating mf poloidal scalar
Base.@propagate_inbounds function s_mf_fac_cs(l, n)
    c1 = (2 * l + 4 * n - 3)
    c2 = -2 * (2 * l + 4 * n - 1)
    c3 = (2 * l + 4 * n + 1)
    fac = 1 / (sqrt(2l * (1 + l) * (-3 + 2 * l + 4 * n) * (-1 + 2 * l + 4 * n) * (1 + 2 * l + 4 * n)))
    return c1, c2, c3, fac
end

Base.@propagate_inbounds function s_mf_pre(js, rls, l, m, n,rgrid::Vector{T}, i) where {T}
    c1, c2, c3, fac = s_mf_fac_cs(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    J = js[l][1, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    Jnm2 = (n - 2 < 0) ? zero(T) : js[l][1, n-1, i]
    return fac * rl * (c1 * J + c2 * Jnm1 + c3 * Jnm2)
end

Base.@propagate_inbounds function d_s_mf_pre(js, rls, l, m, n,rgrid::Vector{T}, i) where {T}
    c1, c2, c3, fac = s_mf_fac_cs(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    Jnm2 = (n - 2 < 0) ? zero(T) : js[l][1, n-1, i]
    dJnm2 = (n - 2 < 0) ? zero(T) : js[l][2, n-1, i]
    return (c1 * d_rlJ(l, r, rl, J, dJ) + c2 * d_rlJ(l, r, rl, Jnm1, dJnm1) + c3 * d_rlJ(l, r, rl, Jnm2, dJnm2)) * fac
end

Base.@propagate_inbounds function d2_s_mf_pre(js, rls, l, m, n,rgrid::Vector{T}, i) where {T}
    c1, c2, c3, fac = s_mf_fac_cs(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    Jnm2 = (n - 2 < 0) ? zero(T) : js[l][1, n-1, i]
    dJnm2 = (n - 2 < 0) ? zero(T) : js[l][2, n-1, i]
    d2Jnm2 = (n - 2 < 0) ? zero(T) : js[l][3, n-1, i]
    return (c1 * d2_rlJ(l, r, rl, J, dJ, d2J) + c2 * d2_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1) + c3 * d2_rlJ(l, r, rl, Jnm2, dJnm2, d2Jnm2)) * fac
end

Base.@propagate_inbounds function d3_s_mf_pre(js, rls, l, m, n,rgrid::Vector{T}, i) where {T}
    c1, c2, c3, fac = s_mf_fac_cs(T(l), T(n))
    r = rgrid[i]
	rl = rls[l][i]
    J = js[l][1, n+1, i]
    dJ = js[l][2, n+1, i]
    d2J = js[l][3, n+1, i]
    d3J = js[l][4, n+1, i]
    Jnm1 = (n - 1 < 0) ? zero(T) : js[l][1, n, i]
    dJnm1 = (n - 1 < 0) ? zero(T) : js[l][2, n, i]
    d2Jnm1 = (n - 1 < 0) ? zero(T) : js[l][3, n, i]
    d3Jnm1 = (n - 1 < 0) ? zero(T) : js[l][4, n, i]
    Jnm2 = (n - 2 < 0) ? zero(T) : js[l][1, n-1, i]
    dJnm2 = (n - 2 < 0) ? zero(T) : js[l][2, n-1, i]
    d2Jnm2 = (n - 2 < 0) ? zero(T) : js[l][3, n-1, i]
    d3Jnm2 = (n - 2 < 0) ? zero(T) : js[l][4, n-1, i]
    return (c1 * d3_rlJ(l, r, rl, J, dJ, d2J, d3J) + c2 * d3_rlJ(l, r, rl, Jnm1, dJnm1, d2Jnm1, d3Jnm1) + c3 * d3_rlJ(l, r, rl, Jnm2, dJnm2, d2Jnm2, d3Jnm2)) * fac
end