module Utils

export appendit!

function appendit!(is, js, aijs, i, j, aij; thresh=sqrt(eps()))
    if !isnothing(aij) && (abs(aij) > thresh)
        push!(is, i)
        push!(js, j)
        push!(aijs, aij)
    end
end
function appendit!(is::Vector{Int}, js::Vector{Int}, aijs::Union{Vector{T},Vector{Complex{T}}}, i::Int, j::Int, aij::T; thresh=sqrt(eps())) where T
    if (abs(aij) > thresh)
        push!(is, i)
        push!(js, j)
        push!(aijs, aij)
    end
end


end #module