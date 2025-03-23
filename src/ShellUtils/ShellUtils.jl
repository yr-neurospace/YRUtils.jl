module ShellUtils

include("../BaseUtils/BaseUtils.jl")
using .BaseUtils

# shell_utils.jl
export para_cmds, find_cmd, pipe_stderr_one, cmd_valid
# file_compress.jl
export recur_pigz

include("base_utils.jl")
include("file_compress.jl")

end # module ShellUtils