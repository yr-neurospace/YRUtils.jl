module ShellUtils

include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

# shell_utils.jl
export para_cmds, find_cmd, pipe_stderr_one, cmd_valid
# svg_compress.jl
export svgcleaner
# file_compress.jl
export recur_pigz
# file_integrity_check.jl
export md5_gen, md5_check

include("base_utils.jl")
include("svg_compress.jl")
include("file_compress.jl")
include("file_integrity_check.jl")

end # module ShellUtils