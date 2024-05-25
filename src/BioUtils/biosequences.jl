const NUCLEOTIDEs = ("A", "C", "G", "T", "U", "N")

abstract type Align end

"""
    check_nucleotide(nucleotide::AbstractString) -> Bool

Check whether `nucleotide` contains valid nucleotide characters: ["A", "C", "G", "T", "U", "N"].
"""
function check_nucleotide(nucleotide::AbstractString)
    if all(in.(unique(convert(Vector{String}, split(nucleotide, ""))), Ref(NUCLEOTIDEs)))
        true
    else
        false
    end
end

"""
    show(aln::Vector{String}, ::Type{Align}; colors::Vector{Tuple{Int64,Int64,Int64}}=Tuple{Int64,Int64,Int64}[]) -> nothing

Show colorful alignments `aln` with `colors` in the format: `[(199, 33, 221), (209, 74, 0), (0, 140, 0), (0, 127, 177)]`.

`(R, G, B)` with each ranging from 0 to 255.
"""
function show(aln::Vector{String}, ::Type{Align}; colors::Vector{Tuple{Int64,Int64,Int64}}=Tuple{Int64,Int64,Int64}[])
    aln = strip.(aln)

    if length(unique(length.(aln))) != 1 || unique(length.(aln)) == 0
        @error "all sequences in aln should have the same length and are non-empty"
    end

    unique_chars = unique(convert(Vector{String}, split(join(aln), "")))
    if length(colors) != length(unique_chars)
        colors = get_distinguishable_colors(length(unique_chars), [RGB(1, 1, 1), RGB(0, 0, 0)]; dropseed=true, encoding=:RGB)
    end

    char_colors = Dict(zip(sort(unique_chars), join.(colors, ";")))
    aln_colored = similar(aln)
    for i in eachindex(aln)
        seq_chars = convert(Vector{String}, split(aln[i], ""))
        seq_colored = join(string.("\033[38;2;", [char_colors[k] for k in seq_chars], "m", seq_chars, "\033[0m"))
        aln_colored[i] = seq_colored
    end

    print(join(aln_colored, "\n"))
end