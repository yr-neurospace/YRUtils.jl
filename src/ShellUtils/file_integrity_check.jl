"""
    md5_gen(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, md5_path::AbstractString="", kwargs...) -> Vector{String}

Generate MD5 sums for files matching given `pattern` in the directory `dir` in a recursive mode if `recursive = true`.

`md5sum` will be found automatically using `which md5sum` if `md5_path` not given.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function md5_gen(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, md5_path::AbstractString="", kwargs...)
    if isempty(dir)
        @error "dir is empty"
    else
        dir = abspath(expanduser(dir))
    end

    if isempty(md5_path)
        md5_path = find_cmd("md5sum")
    end

    cmd_valid(Cmd(string.([md5_path, "--version"])); return_false=false)

    all_files = list_files(dir, pattern; recursive=recursive, full_name=true)

    # for storing the MD5 sum of each file
    md5_res_dict = Dict(k => Dict{String,String}() for k in unique(dirname.(all_files)))
    for (k, v) in zip(dirname.(all_files), all_files)
        md5_res_dict[k][v] = ""
    end

    para_cmds(all_files; kwargs...) do x
        cmd = Cmd(string.([md5_path; x]))
        md5_res_dict[dirname(x)][x] = read(cmd, String)
    end

    for (k, v) in md5_res_dict
        open(joinpath(k, "MD5.txt"), "w") do io
            print(io, join(values(v), ""))
        end
    end

    joinpath.(unique(dirname.(all_files)), "MD5.txt")
end

"""
    md5_check(dir::AbstractString, output_file::AbstractString="", pattern::Regex=r""; recursive::Bool=true, md5_path::AbstractString="", kwargs...) -> String

Check MD5 sums contained in MD5 files matching `pattern` found in the directory `dir` and output the results to `output_file`.

If `output_file` is empty, then the results will be output to `stdout`.

`md5sum` will be found automatically using `which md5sum` if `md5_path` not given.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function md5_check(dir::AbstractString, output_file::AbstractString="", pattern::Regex=r""; recursive::Bool=true, md5_path::AbstractString="", kwargs...)
    if isempty(dir)
        @error "dir is empty"
    else
        dir = abspath(expanduser(dir))
    end

    if isempty(md5_path)
        md5_path = find_cmd("md5sum")
    end

    cmd_valid(Cmd(string.([md5_path, "--version"])); return_false=false)

    if pattern == r""
        pattern = r"MD5\.txt$"
    end

    if isempty(output_file)
        output_io = stdout
    else
        output_file = abspath(expanduser(output_file))
        output_io = open(output_file, "w")
    end

    all_files = list_files(dir, pattern; recursive=recursive, full_name=true)
    lk = ReentrantLock()
    para_cmds(all_files; kwargs...) do x
        # read MD5 sums and corresponding files
        # get full path of each file
        md5_sums = open(x, "r") do io
            ms = match.(r"(\S+) {2}(.+)", readlines(io))
            if isnothing(ms)
                @error string("parse ", x, " failed")
            end
            string.(convert(Vector{String}, getindex.(ms, 1)), "  ",
                joinpath.(dirname(x),
                    convert(Vector{String}, getindex.(ms, 2))))
        end

        # check MD5 sums
        # output all stdout/stderr to a single file
        cmd = Cmd(string.([md5_path, "-w", "-c"]))
        check_res_stdout_io_bf = IOBuffer()
        check_res_stderr_io_bf = IOBuffer()
        try
            run(pipeline(cmd;
                    stdin=IOBuffer(join(md5_sums, "\n")),
                    stdout=check_res_stdout_io_bf,
                    stderr=check_res_stderr_io_bf); wait=true)
        catch e
            @info "in most cases, the following ProcessFailedException captured and printed as a warning only indicates that some files check failed"
            @warn e
        end
        check_res = string(String(take!(check_res_stdout_io_bf)),
            String(take!(check_res_stderr_io_bf)))

        lock(lk) do
            print(output_io, string("==> ", x, ":\n\n", check_res, "\n\n"))
        end
    end

    close(output_io)

    return output_file
end
