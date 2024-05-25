module BioUtils

using ColorTypes
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
export fastqc, trimgalore

include("biosequences.jl")
include("parse_cigar.jl")
include("ngs.jl")

end # module BioUtils