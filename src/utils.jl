module Utils

export appendit!

function appendit!(is, js, aijs, i, j, aij; thresh=sqrt(eps()))
    if !isnothing(aij) && (abs(aij) > thresh)
        push!(is, i)
        push!(js, j)
        push!(aijs, aij)
    end
end



end #module