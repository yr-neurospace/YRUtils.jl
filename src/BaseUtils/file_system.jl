"""
    list_files(dir::AbstractString=pwd(), pattern::Regex=r".+"; recursive::Bool=true, full_name::Bool=true) -> Vector{String}

Return the file names matching `pattern` in the directory `dir` or the current directory if not given in a recursive mode if `recursive = true`.

Set `full_name = true` to get full paths.
"""
function list_files(dir::AbstractString=pwd(), pattern::Regex=r".+"; recursive::Bool=true, full_name::Bool=true)
    if isempty(dir)
        @error "dir is empty"
    else
        abs_dir = abspath(expanduser(dir))
    end

    if pattern == r""
        @error "pattern is empty"
    end

    all_files = String[]
    if recursive
        @info string("search the directory ", abs_dir, " in recursive mode")
        for (root, _, files) in walkdir(abs_dir)
            all_files = [all_files; joinpath.(root, files)]
        end
    else
        for file in readdir(abs_dir; join=true)
            if isfile(file)
                all_files = [all_files; file]
            end
        end
    end
    all_files = all_files[map(x -> !isnothing(x), match.(pattern, all_files))]
    @info string("find ", length(all_files), " files in total")

    return if full_name
        unique(all_files)
    else
        unique(basename.(all_files))
    end
end