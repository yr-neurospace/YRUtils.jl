module BioUtils

using DataFrames

include("../ShellUtils/ShellUtils.jl")
using .ShellUtils
include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

# ngs.jl
export fastqc, trimgalore, auto_detect_fastq_read_type
# ncbi.jl
export ncbi_sra_prefetch, ncbi_sra_dump

include("ngs.jl")
include("ncbi.jl")

end # module BioUtils
