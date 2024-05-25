# valid CIGAR operations (characters)
const CIGAR_OPs = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
const CIGAR_CHARs = ("M", "I", "D", "N", "S", "H", "P", "=", "X", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
const CIGAR_Pattern = r"([0-9]{1,})([MIDNSHP=X]{1})"

"""
    parse_cigar_only(cigar::AbstractString) -> Tuple{String, String}

Parse the CIGAR string `cigar` into aligned pseudo-sequences.
"""
function parse_cigar_only(cigar::AbstractString)
    cigar = strip(cigar)

    if isempty(cigar)
        @error "the argument cigar is empty"
    end
    if !all(in.(unique(convert(Vector{String}, split(cigar, ""))), Ref(CIGAR_CHARs)))
        @error string("CIGAR strings are only allowed to contain the following characters: ", CIGAR_CHARs)
    end

    cigar_match = collect(eachmatch(CIGAR_Pattern, cigar))
    if isempty(cigar_match)
        @error string("parse failed on the CIGAR string: ", cigar)
    end

    ref_align = query_align = ""
    for m in cigar_match
        sub_ref_align, sub_query_align = parse_cigar_op(parse(Int, m[1]), m[2])
        ref_align = join([ref_align, sub_ref_align])
        query_align = join([query_align, sub_query_align])
    end

    if length(ref_align) == length(query_align)
        return (ref_align, query_align)
    else
        @error "the lengthes of reference and query aligned are unequal"
    end
end

"""
    parse_cigar(cigar::AbstractString, ref::AbstractString, query::AbstractString, ref_pos::Int) -> Tuple{String, STring}

Parse the CIGAR string `cigar` into aligned sequences, specified in `ref` (reference) and `query` (query).

`ref_pos` represents **1-based leftmost mapping position of the first CIGAR operation that "consumes" a reference base**,
so `ref_pos` must be greater than or equal to 1 and less than or equal to the length of `ref`.

For more detailed explanations about CIGAR string, see https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf.
"""
function parse_cigar(cigar::AbstractString, ref::AbstractString, query::AbstractString, ref_pos::Int)
    cigar, ref, query = strip(cigar), strip(ref), strip(query)

    if isempty(cigar) || isempty(ref) || isempty(query) || ref_pos < 1
        @error "some of the four arguments cigar, ref, query, and ref_pos is/are invalid (empty or <1)"
    end
    if !all(in.(unique(convert(Vector{String}, split(cigar, ""))), Ref(CIGAR_CHARs)))
        @error string("CIGAR strings are only allowed to contain the following characters: ", CIGAR_CHARs)
    end
    if ref_pos < 1 || ref_pos > length(ref)
        @error string("ref_pos (", ref_pos, ") is either < 1 or > length(ref) (", length(ref), ")")
    end

    cigar_match = collect(eachmatch(CIGAR_Pattern, cigar))
    if isempty(cigar_match)
        @error string("parse failed on the CIGAR string: ", cigar)
    end

    ref_align = query_align = ""
    ref_pointer, query_pointer = ref_pos, 1
    for m in cigar_match
        sub_ref_align, sub_query_align, ref_pointer, query_pointer = parse_cigar_op(parse(Int, m[1]), m[2], ref, query, ref_pointer, query_pointer)
        ref_align = join([ref_align, sub_ref_align])
        query_align = join([query_align, sub_query_align])
    end

    if length(ref_align) == length(query_align)
        return (ref_align, query_align)
    else
        @error "the lengthes of reference and query are unequal"
    end
end

"""
    parse_cigar_op(n::Int, op::AbstractString) -> Tuple{String, STring}

Parse each CIGAR operation into aligned pseudo-sequences.

e.g. for `6M`, `n` = 6, and `op` = "M".
"""
function parse_cigar_op(n::Int, op::AbstractString)
    if op in ("M", "=", "X")  # (consuming query/reference: yes/yes)
        ref = query = join(repeat([op], n))
        return (ref, query)
    elseif op == "I"  # (yes/no)
        ref = join(repeat(["-"], n))
        query = join(repeat([op], n))
        return (ref, query)
    elseif op in ("D", "N")  # (no/yes)
        ref = join(repeat([op], n))
        query = join(repeat(["-"], n))
        return (ref, query)
    elseif op in ("S", "H")  # (yes/no), (no/no)
        ref = join(repeat(["."], n))
        query = join(repeat([op], n))
        return (ref, query)
    elseif op == "P"  # (no/no)
        ref = join(repeat(["-"], n))
        query = join(repeat(["-"], n))
        return (ref, query)
    else
        @error string("unrecognized CIGAR operation: ", op)
    end
end

"""
    parse_cigar_op(n::Int, op::AbstractString, ref::AbstractString, query::AbstractString, ref_pointer::Int, query_pointer::Int) -> Tuple{String, STring}

Parse each CIGAR operation into aligned sequences.

e.g. for `6M`, `n` = 6, and `op` = "M".

`ref` and `query` indicate the reference and query respectively.

`ref_pointer` and `query_pointer` indicate the next indeices consuming the reference and query respectively.
"""
function parse_cigar_op(n::Int, op::AbstractString, ref::AbstractString, query::AbstractString, ref_pointer::Int, query_pointer::Int)
    if op in ("M", "=", "X")  # (consuming query/reference: yes/yes)
        return (ref[ref_pointer:(ref_pointer+n-1)], query[query_pointer:(query_pointer+n-1)],
            ref_pointer + n, query_pointer + n)
    elseif op == "I"  # (yes/no)
        return (join(repeat(["-"], n)), query[query_pointer:(query_pointer+n-1)],
            ref_pointer, query_pointer + n)
    elseif op in ("D", "N")  # (no/yes)
        return (ref[ref_pointer:(ref_pointer+n-1)], join(repeat(["-"], n)),
            ref_pointer + n, query_pointer)
    elseif op == "S"  # (yes/no)
        return (join(repeat(["."], n)), query[query_pointer:(query_pointer+n-1)],
            ref_pointer, query_pointer + n)
    elseif op == "H"  # (no/no)
        return (join(repeat(["."], n)), join(repeat(["."], n)),
            ref_pointer, query_pointer)
    elseif op == "P"  # (no/no)
        return (join(repeat(["-"], n)), join(repeat(["-"], n)),
            ref_pointer, query_pointer)
    else
        @error string("unrecognized CIGAR operation: ", op)
    end
end
