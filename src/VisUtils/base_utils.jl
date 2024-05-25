"""
    get_distinguishable_colors(n::Int, seed::AbstractVector{T}=RGB{N0f8}[]; dropseed::Bool=false, encoding::Symbol=:RGB, kwargs...) where {T<:Color} -> Vector{Tuple{Int64, Int64, Int64}}

Get distinguishable colors with returned value of the form: `[(199, 33, 221), (209, 74, 0), (0, 140, 0), (0, 127, 177)]`.

`(R, G, B)` with each ranging from 0 to 255.

# Examples
```julia-repl
julia> get_distinguishable_colors(6, [RGB(1, 1, 1), RGB(0, 0, 0)]; dropseed=true, encoding=:RGB)
6-element Vector{Tuple{Int64, Int64, Int64}}:
 (199, 33, 221)
 (209, 74, 0)
 (0, 140, 0)
 (0, 127, 177)
 (209, 172, 0)
 (135, 0, 54)
```
"""
function get_distinguishable_colors(n::Int, seed::AbstractVector{T}=RGB{N0f8}[]; dropseed::Bool=false, encoding::Symbol=:RGB, kwargs...) where {T<:Color}
    colors = distinguishable_colors(n, seed; dropseed=dropseed, kwargs...)
    if encoding == :RGB
        map(x -> Int64.(round.((x.r, x.g, x.b) .* 255; digits=0)), convert.(RGB{Float64}, colors))
    else
        @error "only :RGB is supported yet for encoding"
    end
end