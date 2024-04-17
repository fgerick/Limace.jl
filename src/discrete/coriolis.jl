"""
$(TYPEDSIGNATURES)

Equation (113) in Ivers & Phillips (2008).
"""
function C(l,m)
	return (l^2 - 1)*√((l^2 - m^2) / (4l^2 - 1))
end


function _coriolis_TT(lmna, lmnb, r,wr, Ta,Tb; Ω = 2.0)
	l,m,n = lmna
	aij = _inertial_TT(lmna, lmnb, r,wr, Ta,Tb)
    return im*m*Ω/p(l)*aij
end

function _coriolis_SS(lmna, lmnb, r,wr, Sa,Sb; Ω = 2.0)
	l,m,n = lmna
	aij = _inertial_SS(lmna, lmnb, r,wr, Sa,Sb)
    return m*Ω/p(l)*aij
end