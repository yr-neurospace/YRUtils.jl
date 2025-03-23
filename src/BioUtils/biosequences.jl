const DNA_BASE_PAIRS = Dict(
    "A" => "T", "T" => "A",
    "C" => "G", "G" => "C",
    "N" => "N"
)
const RNA_BASE_PAIRS = Dict(
    "A" => "U", "U" => "A",
    "C" => "G", "G" => "C",
    "N" => "N"
)

"""
    rev_seq(seq::AbstractString)

Get reverse sequence of `seq`.
"""
function rev_seq(seq::AbstractString)
    return reverse(seq)
end

"""
    com_dna_seq(seq::AbstractString)

Get complementary sequence of DNA `seq`.
"""
function com_dna_seq(seq::AbstractString)
    return join([DNA_BASE_PAIRS[string(k)] for k in seq])
end

"""
    com_rna_seq(seq::AbstractString)

Get complementary sequence of RNA `seq`.
"""
function com_rna_seq(seq::AbstractString)
    return join([RNA_BASE_PAIRS[string(k)] for k in seq])
end

"""
    rev_com_dna_seq(seq::AbstractString)

Get revrese and complementary sequence of DNA `seq`.
"""
function rev_com_dna_seq(seq::AbstractString)
    return com_dna_seq(rev_seq(seq))
end

"""
    rev_com_rna_seq(seq::AbstractString)

Get revrese and complementary sequence of RNA `seq`.
"""
function rev_com_rna_seq(seq::AbstractString)
    return com_rna_seq(rev_seq(seq))
end
