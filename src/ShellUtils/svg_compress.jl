"""
    svgcleaner(dir::AbstractString; recursive::Bool=true, svgcleaner_path::AbstractString="", kwargs...) -> Vector{String}

Compress SVG files using the **Linux** `svgcleaner` command-line tool (https://github.com/RazrFalcon/svgcleaner).

# Arguments:
- `dir`: the root directory containing SVG files.
- `recursive`: whether to find SVG files recursively.
- `svgcleaner_path`: the path to `svgcleaner`. If not given, it will try to find it using `which svgcleaner`.
- `kwargs`: all other keyword arguments will be passed to `para_cmds` via `kwargs`.

Files found will be returned after processing.

**Note:** make sure that you have `svgcleaner` installed on your computer.

For convenience, you can get the Linux `svgcleaner` package from the `data` directory in package `YRUtils`.
"""
function svgcleaner(dir::AbstractString; recursive::Bool=true, svgcleaner_path::AbstractString="", kwargs...)
    if isempty(dir)
        @error "dir is empty"
    else
        dir = abspath(expanduser(dir))
    end

    if isempty(svgcleaner_path)
        svgcleaner_path = find_cmd("svgcleaner")
    end

    cmd_valid(Cmd(string.([svgcleaner_path, "--version"])); return_false=false)

    all_files = list_files(dir, r".+\.svg$"; recursive=recursive, full_name=true)
    para_cmds(x -> run(Cmd(string.([svgcleaner_path, x, x])), all_files; wait=true); kwargs...)

    all_files
end
