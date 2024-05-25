"""
    recur_pigz(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, pigz_options::AbstractString="-k", pigz_path::AbstractString="", kwargs...) -> Vector{String}

Call **Linux** compression tool `pigz` in julia over each file mattched by `pattern` in the directory `dir` in a recursive mode if `recursive = true`.

The path to `pigz` will be found using `which pigz` if `pigz_path = ""`.

Other `pigz` options can be given by `pigz_options`.

Files found will be returned after processing.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function recur_pigz(dir::AbstractString, pattern::Regex=r".+"; recursive::Bool=true, pigz_options::AbstractString="-k", pigz_path::AbstractString="", kwargs...)
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