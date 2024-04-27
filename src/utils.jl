module Utils

export appendit!

# function appendit!(is, js, aijs, i, j, aij; thresh=sqrt(eps()))
#     if !isnothing(aij) && (abs(aij) > thresh)
#         push!(is, i)
#         push!(js, j)
#         push!(aijs, aij)
#     end
# end
function appendit!(is::Array{Int}, js::Array{Int}, aijs::Array{T}, i::Int, j::Int, aij::T; thresh=sqrt(eps())) where T
    if (abs(aij) > thresh)
        push!(is, i)
        push!(js, j)
        push!(aijs, aij)
    end
end


end #module