"""
    para_cmds(f::Function, arr::AbstractArray, f_args::Union{NamedTuple,Nothing}=nothing; num_jobs::Int=1, poll_secs::Real=1, verbose::Bool=true)::Nothing

Apply function `f` to each element of `arr` concurrently. This function is mainly designed to launch external commands contained in `f` concurrently.

To achieve this goal, you must set `wait = false` when calling the `run` function.

You can provide extra arguments to `f` by passing a named tuple to `f_args`.

The first argument of `f` must accept a single element of `arr`.

You can set the number of jobs launched at the same time via `num_jobs`.

`poll_secs` control the time inteval to poll tasks.

# Examples
```julia-repl
julia> para_cmds(x -> run(Cmd(["sleep", string(x)]); wait=true), [20, 25, 30, 35, 40]; num_jobs=3, poll_secs=3)
[ Info: all jobs done!
```
"""
function para_cmds(f::Function, arr::AbstractArray, f_args::Union{NamedTuple,Nothing}=nothing; num_jobs::Int=1, poll_secs::Real=1, verbose::Bool=true)::Nothing
    arr = copy(arr)

    if isempty(arr)
        throw(ArgumentError("arr is empty"))
    end
    if num_jobs < 1
        @warn "invalid num_jobs, reset to 1"
    end
    if poll_secs < 0
        @warn "invalid poll_secs, reset to 1"
    end
    if length(arr) < num_jobs
        num_jobs = length(arr)
    end

    # initialize tasks
    task_pool = Array{Task}(undef, num_jobs)
    for i in 1:num_jobs
        if !isempty(arr)
            ele = popfirst!(arr)
            if isnothing(f_args)
                task_pool[i] = @async f(ele)
            else
                task_pool[i] = @async f(ele, f_args...)
            end
        else
            break
        end
    end

    # poll and launch new tasks
    while true
        if isempty(arr) && all(istaskdone.(task_pool))
            if verbose
                @info "all jobs done!"
            end
            break
        end

        sleep(poll_secs)

        for i in 1:num_jobs
            if istaskdone(task_pool[i]) && !isempty(arr)
                ele = popfirst!(arr)
                if isnothing(f_args)
                    task_pool[i] = @async f(ele)
                else
                    task_pool[i] = @async f(ele, f_args...)
                end
            end
        end
    end
end


"""
    find_cmd(cmd_name::AbstractString; return_nothing::Bool=false, verbose::Bool=true)::Union{Nothing, AbstractString}

Search the command `cmd_name` path using `which cmd_name`.

If `return_nothing` is `true`, then return `nothing` when `cmd_name` cannot be found instead of throwing an error.
"""
function find_cmd(cmd_name::AbstractString; return_nothing::Bool=false, verbose::Bool=true)::Union{Nothing, AbstractString}
    cmd_path = ""
    try
        cmd_path = open(Cmd(string.(["which", cmd_name])), "r", stdin) do io
            readchomp(io)
        end
    catch e
        if return_nothing
            nothing
        elseif isa(e, ProcessFailedException)
            @error string("command ", cmd_name, " NOT found in the current PATH")
        else
            @warn "an unexpected exception occurred"
            throw(e)
        end
    else
        if verbose
            @info string("command ", cmd_name, " found in path: ", cmd_path)
            cmd_path
        else
            cmd_path
        end
    end
end

"""
    pipe_stderr_one(cmds::Vector{String}, stderr_file::AbstractString; final_stdout_file::AbstractString="")::Expr

Pipe commands and redirect all `stderr`s to a single log file.

optionally, you can redirect the final `stdout` to `final_stdout_file`.

# Examples
```julia-repl
julia> ex = pipe_stderr_one(["cat test.txt", "cut -f 1,3,5", "sort -k1,1 -k2,2"], "log.err", final_stdout_file="result.txt")
quote
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:163 =#
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:163 =# @info "running commands ..."
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:164 =#
    run(pipeline(pipeline(pipeline(pipeline(pipeline(pipeline(`cat test.txt`, stdout=`cut -f 1,3,5`), stderr>Base.FileRedirect("log.err.tmp.cmd1", false)), stdout=`sort -k1,1 -k2,2`), stderr>Base.FileRedirect("log.err.tmp.cmd2", false)), stdout>Base.FileRedirect("result.txt", false)), stderr>Base.FileRedirect("log.err.tmp.cmd3", false)))
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:165 =#
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:165 =# @info "collecting all standard errors and writing them into a single log file ..."
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:166 =#
    run(pipeline(`cat log.err.tmp.cmd1 log.err.tmp.cmd2 log.err.tmp.cmd3`, stdout>Base.FileRedirect("log.err", false)))
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:167 =#
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:167 =# @info "deleting all temporary files ..."
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:168 =#
    run(`rm -rf log.err.tmp.cmd1 log.err.tmp.cmd2 log.err.tmp.cmd3`)
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:169 =#
    #= /home/yangrui/mywd/proj/YRUtils.jl/src/shell_utils.jl:169 =# @info "all jobs done!"
end

julia> eval(ex)
```
"""
function pipe_stderr_one(cmds::Vector{String}, stderr_file::AbstractString; final_stdout_file::AbstractString="")::Expr
    cmd_num = length(cmds)
    if cmd_num < 2
        @error "only support to pipe at least 2 commands"
    end

    cmd_pipes = Vector{Base.AbstractCmd}(undef, cmd_num)
    for i in 1:cmd_num
        if i == 1
            cmd_pipes[i] = pipeline(Cmd(Base.shell_split(cmds[i]));
                stdout=Cmd(Base.shell_split(cmds[i+1])),
                stderr=string(stderr_file, ".tmp.cmd", i))
        elseif i == cmd_num
            cmd_pipes[i] = pipeline(cmd_pipes[i-1];
                stdout=if final_stdout_file != ""
                    final_stdout_file
                else
                    stdout
                end,
                stderr=string(stderr_file, ".tmp.cmd", i))
        else
            cmd_pipes[i] = pipeline(cmd_pipes[i-1];
                stdout=Cmd(Base.shell_split(cmds[i+1])),
                stderr=string(stderr_file, ".tmp.cmd", i))
        end
    end

    cat_tmp_cmd_pipe = pipeline(`cat $(stderr_file).tmp.cmd$(1:cmd_num)`, stdout=stderr_file)

    rm_tmp_cmd = `rm -rf $(stderr_file).tmp.cmd$(1:cmd_num)`

    return quote
        @info ("running commands ...")
        run($(cmd_pipes[cmd_num]))
        @info ("collecting all standard errors and writing them into a single log file ...")
        run($cat_tmp_cmd_pipe)
        @info ("deleting all temporary files ...")
        run($rm_tmp_cmd)
        @info "all jobs done!"
    end
end

"""
    cmd_valid(cmd::Base.AbstractCmd; return_false::Bool=false) -> Bool

Check whether the `cmd` is valid. If not return `false` when `return_false = true` or raise an error.

# Examples
```julia-repl
julia> cmd_valid(`pigz --version`)
pigz 2.6
true

julia> cmd_valid(`pig --version`; return_false=true)
false
```
"""
function cmd_valid(cmd::Base.AbstractCmd; return_false::Bool=false)
    try
        run(cmd)
    catch e
        if isa(e, Base.IOError)
            if return_false
                return false
            else
                @error string("cmd ", cmd, " is NOT valid")
            end
        else
            @warn "an unexpected exception occurred"
            throw(e)
        end
    else
        @info string("cmd ", cmd, " is valid")
        return true
    end
end
