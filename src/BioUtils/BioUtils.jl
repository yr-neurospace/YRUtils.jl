module BioUtils

using DataFrames
using CSV
using Colors
using Humanize

include("../ShellUtils/ShellUtils.jl")
using .ShellUtils
include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

# ngs.jl
export fastqc, trimgalore, auto_detect_fastq_read_type, bam_index, fa_num, fq_num, bam_num, subsample_fastq
# ncbi.jl
export ncbi_sra_prefetch, ncbi_sra_dump
# parse_cigar.jl
export valid_cigar, parse_cigar_op, parse_cigar, show_align
# biosequences.jl
export rev_seq, com_dna_seq, com_rna_seq, rev_com_dna_seq, rev_com_rna_seq

include("ngs.jl")
include("ncbi.jl")
include("parse_cigar.jl")
include("biosequences.jl")

end # module BioUtils
