# ∫₀¹ f(r) r² dr
function rquad(nr) 
	# _r, _wr = gaussjacobi(nr,0.0,2.0) 
	_r, _wr = gausslegendre(nr)
	_r = (_r .+ 1)/2
	_wr /= 2
	# _wr /= 8
	return _r, _wr
end


function ∫dr(f::F,r::Vector{T},wr::Vector{T})::Complex{T} where {F,T}
	out = zero(ComplexF64)
	for (r,w) in zip(r,wr)
		out+=f(r)*r^2*w
	end
	return out
end