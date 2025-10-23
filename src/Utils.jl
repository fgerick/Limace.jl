module Utils

export appendit!

@inline function appendit!(is, js, aijs, lck::ReentrantLock, i, j, aij; thresh=sqrt(eps()))
    if !isnothing(aij) && (abs(aij) > thresh)
        lock(lck) do
            push!(is, i)
            push!(js, j)
            push!(aijs, aij)
        end
    end
end
@inline function appendit!(is::Vector{Int}, js::Vector{Int}, aijs::Union{Vector{T},Vector{Complex{T}}}, lck::ReentrantLock, i::Int, j::Int, aij::T; thresh=sqrt(eps())) where T
    if (abs(aij) > thresh)
        lock(lck) do 
            push!(is, i)
            push!(js, j)
            push!(aijs, aij)
        end
    end
end


end #module