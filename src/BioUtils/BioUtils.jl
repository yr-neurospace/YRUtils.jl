module BioUtils

using DataFrames

include("../ShellUtils/ShellUtils.jl")
using .ShellUtils
include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

# ngs.jl
export fastqc, trimgalore, auto_detect_fastq_read_type

include("ngs.jl")

end # module BioUtils
