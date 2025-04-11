"""
$(TYPEDSIGNATURES)

Computes the advection term for a poloidal/toroidal background flow `U0`, a velocity basis `u`. Equivalent to `-lorentz(u,u,U0)`.
"""
advection(u::Tu, U0::BasisElement{T0,TH,T}; threads=false) where {Tu<:Basis,T0<:Basis,TH<:Helmholtz,T} = -_lorentz(Val(threads), u, u, U0)
