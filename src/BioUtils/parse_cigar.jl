# For various specifications of SAM/BAM and related high-throughput sequencing file formats,
# refer to https://github.com/samtools/hts-specs.
# More specifically, for how to interpret CIGAR strings,
# refer to https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf.

const CIGAR_OPs = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
const CIGAR_CHARs = ("M", "I", "D", "N", "S", "H", "P", "=", "X", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
const CIGAR_UNIT_PATTERN = r"([0-9]{1,})([MIDNSHP=X]{1})"
const CIGAR_STRING_PATTERN = r"^([0-9]{1,}H){0,1}([0-9]{1,}S){0,1}([0-9]{1,}[MIDNP=X]{1}){0,}([0-9]{1,}S){0,1}([0-9]{1,}H){0,1}$"

"""
    valid_cigar(cigar::AbstractString; return_false::Bool=false)

Check whether `cigar` is valid.
"""
function valid_cigar(cigar::AbstractString; return_false::Bool=false)
    if isnothing(match(CIGAR_STRING_PATTERN, cigar))
        if return_false
            return false
        else
            @error "invalid CIGAR string: either contains invalid characters or has malformed format"
        end
    end
end

"""
    parse_cigar_op(n::Int, op::AbstractString)

Parse CIGAR operation `<n><op>`.

Return a tuple: `(<pseudo-reference string>, <pseudo-query string>)`.
"""
function parse_cigar_op(n::Int, op::AbstractString)
    if op == "M"
        # Alignment match (can be a sequence match or mismatch)
        # Consume both reference and query
        ref = query = join(repeat([op], n))
    elseif op == "="
        # Sequence match
        # Consume both reference and query
        ref = query = join(repeat([op], n))
    elseif op == "X"
        # Sequence mismatch
        # Consume both reference and query
        ref = query = join(repeat([op], n))
    elseif op == "I"
        # Insertion to the reference
        # Consume query only
        ref = join(repeat(["-"], n))
        query = join(repeat([op], n))
    elseif op == "D"
        # Deletion from the reference
        # Consume reference only
        ref = join(repeat([op], n))
        query = join(repeat(["-"], n))
    elseif op == "N"
        # Skipped region from the reference
        # Consume reference only
        ref = join(repeat([op], n))
        query = join(repeat(["-"], n))
    elseif op == "S"
        # Soft clipping (clipped sequences present in SEQ)
        # Consume query only
        ref = join(repeat(["."], n))
        query = join(repeat([op], n))
    elseif op == "H"
        # Hard clipping (clipped sequences NOT present in SEQ)
        # Consume neither
        ref = join(repeat(["."], n))
        query = join(repeat([op], n))
    elseif op == "P"
        # Padding (silent deletion from padded reference)
        # Consume neither
        ref = join(repeat(["-"], n))
        query = join(repeat(["-"], n))
    else
        @error string("unrecognized CIGAR operation: ", op)
    end
    return (ref, query)
end

"""
    parse_cigar_op(n::Int, op::AbstractString,
    ref_str::AbstractString, query_str::AbstractString,
    ref_pos::Int, query_pos::Int)

Parse CIGAR operation `<n><op>`.

Return a tuple: `(<sub-reference string>, <sub-query string>, <next reference index>, <next query index>)`.

Note: both `ref_pos` and `query_pos` are 1-based coordinates.
"""
function parse_cigar_op(n::Int, op::AbstractString,
    ref_str::AbstractString, query_str::AbstractString,
    ref_pos::Int, query_pos::Int)
    if op == "M"
        # Alignment match (can be a sequence match or mismatch)
        # Consume both reference and query
        return (
            ref_str[ref_pos:(ref_pos+n-1)],
            query_str[query_pos:(query_pos+n-1)],
            ref_pos + n, query_pos + n
        )
    elseif op == "="
        # Sequence match
        # Consume both reference and query
        return (
            ref_str[ref_pos:(ref_pos+n-1)],
            query_str[query_pos:(query_pos+n-1)],
            ref_pos + n, query_pos + n
        )
    elseif op == "X"
        # Sequence mismatch
        # Consume both reference and query
        return (
            ref_str[ref_pos:(ref_pos+n-1)],
            query_str[query_pos:(query_pos+n-1)],
            ref_pos + n, query_pos + n
        )
    elseif op == "I"
        # Insertion to the reference
        # Consume query only
        return (
            join(repeat(["-"], n)),
            query_str[query_pos:(query_pos+n-1)],
            ref_pos, query_pos + n
        )
    elseif op == "D"
        # Deletion from the reference
        # Consume reference only
        return (
            ref_str[ref_pos:(ref_pos+n-1)],
            join(repeat(["-"], n)),
            ref_pos + n, query_pos
        )
    elseif op == "N"
        # Skipped region from the reference
        # Consume reference only
        return (
            ref_str[ref_pos:(ref_pos+n-1)],
            join(repeat(["-"], n)),
            ref_pos + n, query_pos
        )
    elseif op == "S"
        # Soft clipping (clipped sequences present in SEQ)
        # Consume query only
        return (
            join(repeat(["."], n)),
            query_str[query_pos:(query_pos+n-1)],
            ref_pos, query_pos + n
        )
    elseif op == "H"
        # Hard clipping (clipped sequences NOT present in SEQ)
        # Consume neither
        return (
            join(repeat(["."], n)),
            join(repeat(["."], n)),
            ref_pos, query_pos
        )
    elseif op == "P"
        # Padding (silent deletion from padded reference)
        # Consume neither
        return (
            join(repeat(["."], n)),
            join(repeat(["."], n)),
            ref_pos, query_pos
        )
    else
        @error string("unrecognized CIGAR operation: ", op)
    end
end

"""
    parse_cigar(cigar::AbstractString)

Parse CIGAR string `cigar` into aligned pseudo-reference and pseudo-query sequences.

Return a tuple: `(<pseudo-reference sequence>, <pseudo-query sequence>)`.
"""
function parse_cigar(cigar::AbstractString)
    valid_cigar(cigar; return_false=false)

    cigar_matches = collect(eachmatch(CIGAR_UNIT_PATTERN, cigar))
    ref_str = query_str = ""
    for cigar_match in cigar_matches
        sub_ref_str, sub_query_str = parse_cigar_op(parse(Int, cigar_match[1]), cigar_match[2])
        ref_str = join([ref_str, sub_ref_str])
        query_str = join([query_str, sub_query_str])
    end

    if length(ref_str) == length(query_str)
        return (ref_str, query_str)
    else
        @error "the sequnece lengths of reference and query parsed are unequal"
    end
end

"""
    parse_cigar(cigar::AbstractString,
    ref_str::AbstractString, query_str::AbstractString,
    ref_start_pos::Int; truncate_ref::Bool=true)

Parse CIGAR string `cigar` into aligned reference and query sequences.

Return a tuple: `(<reference sequence>, <query sequence>)`.
"""
function parse_cigar(cigar::AbstractString,
    ref_str::AbstractString, query_str::AbstractString,
    ref_start_pos::Int; truncate_ref::Bool=true)
    valid_cigar(cigar; return_false=false)

    if isempty(ref_str) || isempty(query_str) || ref_start_pos < 1 || ref_start_pos > length(ref_str)
        @error "some of the three arguments ref_str, query_str, and ref_start_pos are invalid"
    end

    cigar_matches = collect(eachmatch(CIGAR_UNIT_PATTERN, cigar))
    ref_pos, query_pos = ref_start_pos, 1
    ref_aln = query_aln = ""
    for cigar_match in cigar_matches
        sub_ref_aln, sub_query_aln, ref_pos, query_pos = parse_cigar_op(parse(Int, cigar_match[1]), cigar_match[2],
            ref_str, query_str, ref_pos, query_pos)
        ref_aln = join([ref_aln, sub_ref_aln])
        query_aln = join([query_aln, sub_query_aln])
    end

    if !truncate_ref
        if ref_start_pos > 1
            left_unaligned_ref_str = ref_str[1:ref_start_pos-1]
            ref_aln = join([left_unaligned_ref_str; ref_aln])
            query_aln = join([repeat([" "], length(left_unaligned_ref_str)); query_aln])
        end
        if ref_pos <= length(ref_str)
            right_unaligned_ref_str = ref_str[ref_pos:length(ref_str)]
            ref_aln = join([ref_aln; right_unaligned_ref_str])
            query_aln = join([query_aln; repeat([" "], length(right_unaligned_ref_str))])
        end
    end

    if length(ref_aln) == length(query_aln)
        return (ref_aln, query_aln)
    else
        @error "the sequnece lengths of reference and query parsed are unequal"
    end
end

"""
    show_align(aligns::Dict{String,Vector{String}}, html_file::AbstractString="";
    seq_font_weight::AbstractString="bold",
    title_font_weight::AbstractString="bold",
    seq_font_family::AbstractString="Courier New",
    title_font_family::AbstractString="Arial",
    seq_font_size::AbstractString="14px",
    title_font_size::AbstractString="12px",
    single_color::AbstractString="#FFFFFF",
    mapping_colors::Dict{String,String}=Dict{String,String}(),
    mapping_background::Bool=true,
    wrap_width::Int=60)

Show aligned sequences in HTML format.

`seq_font_family` should be monospaced font family, such as Courier, Courier New, Consolas,  and Lucida Console.

If `html_file` is not given, HTML content will be printed to `stdout`.
"""
function show_align(aligns::Dict{String,Vector{String}}, html_file::AbstractString="";
    seq_font_weight::AbstractString="bold",
    title_font_weight::AbstractString="bold",
    seq_font_family::AbstractString="Courier New",
    title_font_family::AbstractString="Arial",
    seq_font_size::AbstractString="14px",
    title_font_size::AbstractString="12px",
    single_color::AbstractString="#FFFFFF",
    mapping_colors::Dict{String,String}=Dict{String,String}(),
    mapping_background::Bool=true,
    wrap_width::Int=60)
    uniq_chars = unique(flatten_array(split.(flatten_array(collect(values(aligns))), "")))
    if sort(uniq_chars) != sort(collect(keys(mapping_colors)))
        @info "mapping_colors is invalid, and built-in colors are assigned to it"
        mapping_colors = Dict(zip(uniq_chars, string.("#", hex.(distinguishable_colors(length(uniq_chars), [RGB(1, 1, 1), RGB(0, 0, 0)]; dropseed=true)))))
    end

    function span_str(str::AbstractString)
        if mapping_background
            join(string.("<span style='color: ", single_color, "; background-color: ",
                [mapping_colors[c] for c in split(str, "")], ";'>", split(str, ""), "</span>"))
        else
            join(string.("<span style='color: ", [mapping_colors[c] for c in split(str, "")], "; background-color: ",
                single_color, ";'>", split(str, ""), "</span>"))
        end
    end

    box_divs = ""
    for id in keys(aligns)
        ref_strs, query_strs = collect.(Iterators.partition.(aligns[id], wrap_width))
        aln_span = []
        for str_pair in collect(zip(ref_strs, query_strs))
            aln_span = [aln_span; join(span_str.(str_pair), "<br>")]
        end
        aln_pre = string("<pre class='alnseq'>", join(aln_span, "<br><br>"), "</pre>")
        id_div = string("<div class='alnseq-title'>", id, ": </div>")
        box_divs = string(box_divs, string("<div class='alnseq-container'>", id_div, aln_pre, "</div>"))
    end

    html_doc = string(
        "<!DOCTYPE html>",
        "<html lang='en'>",
        "<head>",
        "<meta charset='UTF-8'>",
        "<title>Aligned Sequences</title>",
        "<style>",
        ".alnseq-container { margin: 0px; padding: 10px; }",
        string(
            ".alnseq-title { font-family: ", title_font_family, ", Courier, monospace; ",
            "font-weight: ", title_font_weight, "; ",
            "font-size: ", title_font_size, "; ",
            "margin-bottom: 5px; }"
        ),
        string(
            ".alnseq { border: 2px solid skyblue; ",
            "border-radius: 0px; ",
            "padding: 10px; ",
            "margin: 0px; ",
            "font-family: ", seq_font_family, "; ",
            "font-weight: ", seq_font_weight, "; ",
            "font-size: ", seq_font_size, "; }"
        ),
        "</style>",
        "</head>",
        "<body>",
        box_divs,
        "</body>"
    )

    if !isempty(html_file)
        open(html_file, "w") do io
            println(io, html_doc)
        end
    else
        println(stdout, html_doc)
    end
end
