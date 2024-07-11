module BioUtils

using ColorTypes, CSV, DataFrames, FASTX, CodecZlib

include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

include("../ShellUtils/ShellUtils.jl")
using .ShellUtils

include("../VisUtils/VisUtils.jl")
using .VisUtils

# biosequences.jl
export check_nucleotide, show
# parse_cigar.jl
export parse_cigar, parse_cigar_only, parse_cigar_op
# ngs.jl
export fastqc, trimgalore, auto_detect_fastq_read_type
# 10x.jl
export rename_10x_fastq_name, run_10x_cellranger_arc_count, run_10x_cellranger_count
# ncbi.jl
export ncbi_sra_prefetch, ncbi_sra_dump

include("biosequences.jl")
include("parse_cigar.jl")
include("ngs.jl")
include("10x.jl")
include("ncbi.jl")

end # module BioUtils