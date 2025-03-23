"""
    flatten_array(arr::AbstractArray)::AbstractArray

Flatten a nested array `arr`.
"""
function flatten_array(arr::AbstractArray)::AbstractArray
    mapreduce(x -> isempty(x) || !(typeof(x) <: AbstractArray) ? x : flatten_array(x), vcat, arr, init=[])
end