"""
    recur_pigz(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, pigz_options::AbstractString="-k", pigz_path::AbstractString="", kwargs...)::AbstractArray

Call **Linux** compression tool `pigz` in julia over each file mattched by `pattern` in the directory `dir` in a recursive mode if `recursive = true`.

The path to `pigz` will be found using `which pigz` if `pigz_path = ""`.

Other `pigz` options can be given by `pigz_options`.

Files found will be returned after processing.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function recur_pigz(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, pigz_options::AbstractString="-k", pigz_path::AbstractString="", kwargs...)::AbstractArray
    if isempty(dir)
        @error "dir is empty"
    else
        dir = abspath(expanduser(dir))
    end

    if isempty(pigz_path)
        pigz_path = find_cmd("pigz")
    end

    cmd_valid(Cmd(string.([pigz_path, "--version"])); return_false=false)

    pigz_options = strip(pigz_options)
    cmd_vec = [pigz_path]
    if !isempty(pigz_options)
        cmd_vec = [cmd_vec; Base.shell_split(pigz_options)]
    end

    all_files = list_files(dir, pattern; recursive=recursive, full_name=true)
    para_cmds(all_files; kwargs...) do x
        cmd = Cmd(string.([cmd_vec; x]))
        @info string("running ", cmd, " ...")
        run(cmd; wait=true)
        @info string("running ", cmd, " done!")
    end

    all_files
end

"""
    pigz(from::AbstractString, to::AbstractString=""; n_threads::Int=0, decompress::Bool=false, keep::Bool=true, pigz_path::AbstractString="")

A wrapper of linux command line tool `pigz`.

If `to` is not empty, then compressed/decompressed output will be redirected to `to`.

If `n_threads <= 0`, then use the default setting of `pigz`.
"""
function pigz(from::AbstractString, to::AbstractString=""; n_threads::Int=0, decompress::Bool=false, keep::Bool=true, pigz_path::AbstractString="")
    from = strip(from)
    to = strip(to)
    pigz_path = strip(pigz_path)

    if isempty(pigz_path)
        pigz_path = find_cmd("pigz")
    end
    cmd_valid(Cmd(string.([pigz_path, "--version"])); return_false=false)

    if isempty(from)
        @error "from cannot be empty"
    else
        from = abspath(expanduser(from))
    end

    if !isempty(to)
        to = abspath(expanduser(to))
    end

    cmd = [pigz_path]
    if keep
        cmd = [cmd; "-k"]
    end
    if n_threads > 0
        cmd = [cmd; "-p"; n_threads]
    end

    if decompress
        if isempty(to)
            cmd = Cmd(string.([cmd; "-d"; from]))
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        else
            cmd = pipeline(Cmd(string.([cmd; "-c"; "-d"; from])); stdout=to)
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        end
    else
        if isempty(to)
            cmd = Cmd(string.([cmd; from]))
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        else
            cmd = pipeline(Cmd(string.([cmd; "-c"; from])); stdout=to)
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        end
    end
end

"""
    xz(from::AbstractString, to::AbstractString=""; n_threads::Int=0, decompress::Bool=false, keep::Bool=true, xz_path::AbstractString="")

A wrapper of linux command line tool `xz`.

If `to` is not empty, then compressed/decompressed output will be redirected to `to`.

If `n_threads <= 0`, then use the default setting of `xz`.
"""
function xz(from::AbstractString, to::AbstractString=""; n_threads::Int=0, decompress::Bool=false, keep::Bool=true, xz_path::AbstractString="")
    from = strip(from)
    to = strip(to)
    xz_path = strip(xz_path)

    if isempty(xz_path)
        xz_path = find_cmd("xz")
    end
    cmd_valid(Cmd(string.([xz_path, "--version"])); return_false=false)

    if isempty(from)
        @error "from cannot be empty"
    else
        from = abspath(expanduser(from))
    end

    if !isempty(to)
        to = abspath(expanduser(to))
    end

    cmd = [xz_path]
    if keep
        cmd = [cmd; "-k"]
    end
    if n_threads > 0
        cmd = [cmd; "-T"; n_threads]
    end

    if decompress
        if isempty(to)
            cmd = Cmd(string.([cmd; "-d"; from]))
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        else
            cmd = pipeline(Cmd(string.([cmd; "-c"; "-d"; from])); stdout=to)
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        end
    else
        if isempty(to)
            cmd = Cmd(string.([cmd; from]))
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        else
            cmd = pipeline(Cmd(string.([cmd; "-c"; from])); stdout=to)
            @info string("running ", cmd, " ...")
            run(cmd; wait=true)
        end
    end
end

"""
    create_tar(froms::Vector{AbstractString}, to::AbstractString; tar_path::AbstractString="")

A wrapper of linux command line tool `tar -cf`.
"""
function create_tar(froms::Vector{String}, to::AbstractString; tar_path::AbstractString="")
    froms = strip.(froms)
    froms = froms[.!isempty.(froms)]
    to = strip(to)
    tar_path = strip(tar_path)

    if isempty(tar_path)
        tar_path = find_cmd("tar")
    end
    cmd_valid(Cmd(string.([tar_path, "--version"])); return_false=false)

    if isempty(froms)
        @error "froms cannot be empty"
    end

    if isempty(to)
        @error "to cannot be empty"
    end

    cmd = Cmd(string.([tar_path; "-cf"; to; froms]))
    @info string("running ", cmd, " ...")
    run(cmd; wait=true)
end

"""
    extract_tar(from::AbstractString, to::AbstractString; tar_path::AbstractString="")

A wrapper of linux command line tool `tar -xf`.

`from` should be an archive file, and `to` should be a valid directory.
"""
function extract_tar(from::AbstractString, to::AbstractString; tar_path::AbstractString="")
    from = strip(from)
    to = strip(to)
    tar_path = strip(tar_path)

    if isempty(tar_path)
        tar_path = find_cmd("tar")
    end
    cmd_valid(Cmd(string.([tar_path, "--version"])); return_false=false)

    if isempty(from)
        @error "from cannot be empty"
    end

    if isempty(to)
        @error "to cannot be empty"
    end

    cmd = Cmd(string.([tar_path; "-xf"; from; "-C"; to]))
    @info string("running ", cmd, " ...")
    run(cmd; wait=true)
end
