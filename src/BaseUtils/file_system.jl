"""
    list_files(dir::AbstractString=pwd(), pattern::Regex=r".+"; recursive::Bool=true, full_name::Bool=true)::Vector{String}

Return the file names matching `pattern` in the directory `dir` or the current directory if not given in a recursive mode if `recursive = true`.

Set `full_name = true` to get full paths.
"""
function list_files(dir::AbstractString=pwd(), pattern::Regex=r".+"; recursive::Bool=true, full_name::Bool=true)::Vector{String}
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

"""
    file_transfer(froms::Vector{String}, tos::Vector{String}; ops::Vector{String}=["soft", "hard", "copy", "move"])::Nothing

Transfer files via building soft/hard links, copying files, or moving files.

Each files in `froms` must be existed, and the directories contained in `tos` must also be existed.

This function will try each mode given in `ops` in order in case the former failed.
"""
function file_transfer(froms::Vector{String}, tos::Vector{String}; ops::Vector{String}=["soft", "hard", "copy", "move"])::Nothing
    OPs = ("soft", "hard", "copy", "move")

    froms = abspath.(expanduser.(froms))
    tos = abspath.(expanduser.(tos))

    if isempty(ops)
        @error "ops cannot be empty"
    end
    if !all(in.(ops, Ref(OPs)))
        @error "some operations are invalid in ops"
    end

    if isempty(froms) || isempty(tos)
        @error "both froms and tos cannot be empty"
    end
    if !all(isfile.(froms))
        @error "some files are invalid in froms"
    end
    if !all(isdir.(dirname.(tos)))
        @error "some directories are invalid in tos"
    end
    if length(froms) != length(tos)
        @error "the length of froms is not equal to the length of tos"
    end

    fail_flag = true
    for op in ops
        op_status = repeat([false], length(tos))

        @info string("try running file_tansfer in ", op, " mode ...")
        for (i, vals) in enumerate(zip(froms, tos))
            from, to = vals
            try
                if op == "soft"
                    symlink(from, to; dir_target=false)
                elseif op == "hard"
                    hardlink(from, to)
                elseif op == "copy"
                    cp(from, to; force=false, follow_symlinks=false)
                elseif op == "move"
                    mv(from, to; force=false)
                end
            catch e
                break
            else
                op_status[i] = true
            end
        end

        if all(op_status)
            fail_flag = false
            @info string("running file_tansfer in ", op, " mode successfully")
            break
        else
            fail_flag = true
            rm.(ops[op_status]; force=false, recursive=false)
            @info string("running file_tansfer in ", op, " mode failed")
        end
    end

    if fail_flag
        @error "cannot run file_transfer in any given mode successfully"
    end
end
