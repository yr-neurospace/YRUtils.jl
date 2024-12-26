"""
    md5_gen(files::Vector{String}, io::IO=stdout, root_dir::AbstractString="")::Vector{Tuple{String, String}}

Generate MD5 checksums for `files`.

You can trim the common path of `files` by specifying `root_dir`.
"""
function md5_gen(files::Vector{String}, io::Union{Nothing,IO}=stdout, root_dir::AbstractString="")::Vector{Tuple{String,String}}
    if isempty(files)
        @error "files is empty"
    end
    root_dir = strip(root_dir)

    md5_arr = Vector{Vector{String}}(undef, length(files))
    Threads.@threads for i in 1:length(files)
        md5_arr[i] = [bytes2hex(open(md5, files[i])), files[i]]
    end

    if !isempty(root_dir)
        md5_tup = collect(zip(first.(md5_arr), replace.(last.(md5_arr), Regex(string("^", rstrip(root_dir, '/'), "/")) => "")))
    else
        md5_tup = collect(zip(first.(md5_arr), last.(md5_arr)))
    end

    if !isnothing(io)
        print(io, join(join.(md5_tup, "  "), "\n"))
    end

    md5_tup
end

"""
    md5_gen(files::Vector{String}, md5_file::AbstractString, root_dir::AbstractString="")::Vector{Tuple{String,String}}

Generate MD5 checksums for `files`.

You can save the MD5 checksums to a file, specified by `md5_file`.

You can trim the common path of `files` by specifying `root_dir`.
"""
function md5_gen(files::Vector{String}, md5_file::AbstractString, root_dir::AbstractString="")::Vector{Tuple{String,String}}
    if isempty(strip(md5_file))
        @error "md5_file cannot be empty"
    end

    open(md5_file, "w") do io
        md5_gen(files, io, root_dir)
    end
end

"""
    md5_check(md5_file::AbstractString, io::Union{Nothing,IO}=stdout;
    ignore_missing_files::Bool=false, add_root_dir::Bool=true, quiet::Bool=true)::Vector{Tuple{String,String}}

Check MD5 checksums reading from `md5_file`.

If you want to skip files not existed, set `ignore_missing_files = true`.

If the current working directory is not exactly the same as the one in which MD5 file is, 
please set `add_root_dir = true`. This will add the path `dirname(md5_file)` to each file within MD5 file.

By default, `quiet = true`, which means that it will set the checking status to "NO" instead of raising an error 
though some MD5 checksums checking failed.
"""
function md5_check(md5_file::AbstractString, io::Union{Nothing,IO}=stdout;
    ignore_missing_files::Bool=false, add_root_dir::Bool=true, quiet::Bool=true)::Vector{Tuple{String,String}}
    if isempty(strip(md5_file))
        @error "md5_file cannot be empty"
    end

    md5_arr = open(md5_file, "r") do io
        md5_str = strip.(readlines(io))
        convert(Vector{Vector{String}}, split.(md5_str[.!isempty.(md5_str)], r"[ *]+"))
    end

    if !all(length.(md5_arr) .== 2)
        @error "invalid MD5 file format"
    end

    if ignore_missing_files
        md5_arr = md5_arr[isfile.(last.(md5_arr))]
    end

    if add_root_dir
        md5_dict = Dict(zip(joinpath.(dirname(md5_file), last.(md5_arr)), first.(md5_arr)))
    else
        md5_dict = Dict(zip(last.(md5_arr), first.(md5_arr)))
    end

    md5_arr_new = md5_gen([file for file in keys(md5_dict)], nothing)
    md5_dict_new = Dict(zip(last.(md5_arr_new), first.(md5_arr_new)))

    md5_check_dict = Dict(file => "" for file in keys(md5_dict))
    for file in keys(md5_dict)
        md5_check_dict[file] = if md5_dict[file] == md5_dict_new[file]
            "OK"
        else
            "NO"
        end
    end

    md5_check_tup = [(if add_root_dir
            replace(file, Regex(string("^", dirname(md5_file), "/")) => "")
        else
            file
        end, check_status) for (file, check_status) in md5_check_dict]

    if !quiet && !all(last.(md5_check_tup) .== "OK")
        @error "MD5 checksums checking failed for some files"
    end

    if !isnothing(io)
        print(io, join(join.(md5_check_tup, ": "), "\n"))
    end

    md5_check_tup
end

"""
    md5_check(md5_file::AbstractString, md5_check_file::AbstractString; kwargs...)::Vector{Tuple{String,String}}

Check MD5 checksums reading from `md5_file` and save the results to a file specified in `md5_check_file`.

`kwargs` will be passed to `md5_check` accepting an `IO` object.
"""
function md5_check(md5_file::AbstractString, md5_check_file::AbstractString; kwargs...)::Vector{Tuple{String,String}}
    if isempty(strip(md5_check_file))
        @error "md5_check_file cannot be empty"
    end

    open(md5_check_file, "w") do io
        md5_check(md5_file, io; kwargs...)
    end
end
