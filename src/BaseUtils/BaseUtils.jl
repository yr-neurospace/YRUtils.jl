module BaseUtils

using MD5

# file misc.jl
export flatten_array
# file_system.jl
export list_files, file_transfer
# file_integrity_check.jl
export md5_gen, md5_check

include("misc.jl")
include("file_system.jl")
include("file_integrity_check.jl")

end # module BaseUtils