const FASTQ_PATTERN = r"\.(fastq|fq|fastq\.gz|fq\.gz)$"

"""
    fastqc(dir::AbstractString, outdir::AbstractString=pwd(), pattern::Regex=r""; recursive::Bool=true, nthreads::Int=1, fastqc_path::AbstractString="", fastqc_options::AbstractString="", multiqc_path::AbstractString="", kwargs...) -> Vector{String}

Call `fastqc` in Julia on FASTQ files matching `pattern` by searching the directory `dir` in a recursive mode if `recursive = true` followed by calling `multiqc` on the results of `fastqc`.

The path to `fastqc`/`multiqc` will be found using `which fastqc`/`which multiqc` if `fastqc_path = ""`/`multiqc_path = ""`. 

Other `fastqc` options can be given by `fastqc_options`.

Files found will be returned after processing.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function fastqc(dir::AbstractString, outdir::AbstractString=pwd(), pattern::Regex=r""; recursive::Bool=true, nthreads::Int=1, fastqc_path::AbstractString="", fastqc_options::AbstractString="", multiqc_path::AbstractString="", kwargs...)
    if isempty(dir) || isempty(outdir)
        @error "dir and/or outdir is empty"
    else
        dir = abspath(expanduser(dir))
        outdir = abspath(expanduser(outdir))
    end

    if isempty(fastqc_path)
        fastqc_path = find_cmd("fastqc")
    end

    cmd_valid(Cmd(string.([fastqc_path, "--version"])); return_false=false)

    if isempty(multiqc_path)
        multiqc_path = find_cmd("multiqc")
    end

    cmd_valid(Cmd(string.([multiqc_path, "--version"])); return_false=false)

    if pattern == r""
        pattern = r"\.(fastq|fq|fastq\.gz|fq\.gz)$"
    end

    if nthreads < 1
        @warn "invalid nthreads, reset to 1"
    end

    fastqc_outdir = joinpath(outdir, "fastqc")
    multiqc_outdir = joinpath(outdir, "multiqc")
    mkpath(fastqc_outdir)
    mkpath(multiqc_outdir)

    fastqc_options = strip(fastqc_options)
    cmd_vec = [fastqc_path, "--threads", nthreads, "--outdir", fastqc_outdir]
    if !isempty(fastqc_options)
        cmd_vec = [cmd_vec; Base.shell_split(fastqc_options)]
    end

    all_files = list_files(dir, pattern; recursive=recursive, full_name=true)
    @info "start running fastqc ..."
    para_cmds(all_files; kwargs...) do x
        stdout_io_bf = IOBuffer()
        stderr_io_bf = IOBuffer()
        cmd = Cmd(string.([cmd_vec; x]))
        run(pipeline(cmd; stdout=stdout_io_bf, stderr=stderr_io_bf); wait=true)
        stdout_stderr_str = string(String(take!(stdout_io_bf)),
            String(take!(stderr_io_bf)), "\n\n")
        print(stdout_stderr_str)
    end

    @info "start running multiqc ..."
    run(Cmd(string.([multiqc_path, "-o", multiqc_outdir, fastqc_outdir])); wait=true)

    all_files
end

"""
    trimgalore(dir::AbstractString, outdir::AbstractString=pwd(), r1_suffix::AbstractString="_1.fq.gz", r2_suffix::AbstractString="_2.fq.gz", pattern::Regex=r""; nthreads::Int=1, recursive::Bool=true, trimgalore_path::AbstractString="", trimgalore_options::AbstractString="--paired --phred33 --quality 20 --length 30 --trim-n", kwargs...) -> Vector{String} or Vector{Tuple{String, String}}

Call `trim_galore` in Julia on FASTQ files matching `r1_suffix` and/or `r2_suffix` by searching the directory `dir` in a recursive mode if `recursive = true`.

If the reads are single-end, then just `r1_suffix` needed.

`pattern` is only used to match the common FASTQ files suffix independent of the type of endedness. The default will work well for most cases, so just leave it empty.

The path to `trim_galore` will be found using `which trim_galore` if `trimgalore_path = ""`. 

Other `trim_galore` options can be given by `trimgalore_options`.

Files found will be returned after processing.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function trimgalore(dir::AbstractString, outdir::AbstractString=pwd(), r1_suffix::AbstractString="_1.fq.gz", r2_suffix::AbstractString="_2.fq.gz", pattern::Regex=r""; nthreads::Int=1, recursive::Bool=true, trimgalore_path::AbstractString="", trimgalore_options::AbstractString="--paired --phred33 --quality 20 --length 30 --trim-n", kwargs...)
    if isempty(dir) || isempty(outdir)
        @error "dir and/or outdir is empty"
    else
        dir = abspath(expanduser(dir))
        outdir = abspath(expanduser(outdir))
    end

    if pattern == r""
        pattern = r"\.(fastq|fq|fastq\.gz|fq\.gz)$"
    end

    if isempty(trimgalore_path)
        trimgalore_path = find_cmd("trim_galore")
    end

    cmd_valid(Cmd(string.([trimgalore_path, "--version"])); return_false=false)

    cutadapt_path = find_cmd("cutadapt")
    pigz_path = find_cmd("pigz")
    cmd_valid(Cmd(string.([cutadapt_path, "--version"])); return_false=false)
    cmd_valid(Cmd(string.([pigz_path, "--version"])); return_false=false)

    if nthreads < 1
        @warn "invalid nthreads, reset to 1"
    end

    mkpath(outdir)

    trimgalore_options = strip(trimgalore_options)
    cmd_vec = [trimgalore_path]
    if !isempty(trimgalore_options)
        cmd_vec = [cmd_vec; Base.shell_split(trimgalore_options)]
    end
    cmd_vec = [cmd_vec; "--cores"; nthreads; "--output_dir"; outdir]

    r1_suffix_pattern = Regex(string(replace(r1_suffix, "." => raw"\."), raw"$"))
    # r2_suffix_pattern = Regex(string(replace(r2_suffix, "." => raw"\."), raw"$"))
    all_r1_files = list_files(dir, r1_suffix_pattern; recursive=recursive, full_name=true)
    if "--paired" in cmd_vec
        all_r2_files = replace.(all_r1_files, r1_suffix_pattern => r2_suffix)
        all_r1_r2_files = collect(zip(all_r1_files, all_r2_files))
        @info "start running trim_galore in paired-end mode ..."
        para_cmds(all_r1_r2_files; kwargs...) do x
            stdout_io_bf = IOBuffer()
            stderr_io_bf = IOBuffer()
            cmd = Cmd(string.([cmd_vec; x[1]; x[2]]))
            @info string("running ", cmd, " ...")
            run(pipeline(cmd; stdout=stdout_io_bf, stderr=stderr_io_bf); wait=true)
            stdout_stderr_str = string(String(take!(stdout_io_bf)),
                String(take!(stderr_io_bf)), "\n\n")
            print(stdout_stderr_str)
        end
    else
        @info "start running trim_galore in single-end mode ..."
        para_cmds(all_r1_files; kwargs...) do x
            stdout_io_bf = IOBuffer()
            stderr_io_bf = IOBuffer()
            cmd = Cmd(string.([cmd_vec; x]))
            @info string("running ", cmd, " ...")
            run(pipeline(cmd; stdout=stdout_io_bf, stderr=stderr_io_bf); wait=true)
            stdout_stderr_str = string(String(take!(stdout_io_bf)),
                String(take!(stderr_io_bf)), "\n\n")
            print(stdout_stderr_str)
        end
    end

    @info "rename FASTQ files ..."
    trimmed_fq_files = list_files(outdir, pattern; recursive=false, full_name=true)
    renamed_trimmed_fq_files = replace.(trimmed_fq_files, Regex(string("(_val_[12]|_trimmed)(?=", pattern.pattern, ")")) => "")
    if isempty(renamed_trimmed_fq_files)
        @error "fail to rename FASTQ files"
    end
    for (s, t) in zip(trimmed_fq_files, renamed_trimmed_fq_files)
        @info string("move ", s, " to ", t, " ...")
        mv(s, t; force=false)
    end

    return if "--paired" in cmd_vec
        all_r1_r2_files
    else
        all_r1_files
    end
end

"""
    auto_detect_fastq_read_type(dir::AbstractString; recursive=true) -> Vector

Detect FASTQ read types and categorize each FASTQ file based on its ID, read type, and rep.
"""
function auto_detect_fastq_read_type(dir::AbstractString; recursive=true)
    fq_files = list_files(dir, FASTQ_PATTERN; recursive=recursive, full_name=true)

    fq_ms = match.(Regex(string("(.+)", "(_rep\\d+){1}", "(_part\\d+){0,1}", "(_[12]|\\.R[12]){0,1}", FASTQ_PATTERN.pattern)), basename.(fq_files))

    if any(isnothing.(fq_ms))
        @error "extract fields from some FASTQ file names failed"
    end

    full_fq_ms = convert(Vector{NTuple{6,Union{Nothing,String}}}, [(m[1], m[2], m[3], m[4], m[5], fq) for (m, fq) in zip(fq_ms, fq_files)])

    if any(isnothing.(getindex.(full_fq_ms, 1))) || any(isnothing.(getindex.(full_fq_ms, 2)))
        @error "some FASTQ file names lack standard ID and/or rep fields"
    end

    if all(.!isnothing.(getindex.(full_fq_ms, 4)))
        @info "it seems that all FASTQ files are paired-end well"
        endedness = "paired"
        fqs_dict = Dict(id => Dict(r => Dict{String,Vector{String}}() for r in unique(unique(getindex.(full_fq_ms, 4)))) for id in unique(getindex.(full_fq_ms, 1)))
        for full_fq_m in full_fq_ms
            if !haskey(fqs_dict[full_fq_m[1]][full_fq_m[4]], full_fq_m[2])
                fqs_dict[full_fq_m[1]][full_fq_m[4]][full_fq_m[2]] = [full_fq_m[6]]
            else
                push!(fqs_dict[full_fq_m[1]][full_fq_m[4]][full_fq_m[2]], full_fq_m[6])
            end
        end
    elseif all(isnothing.(getindex.(full_fq_ms, 4)))
        @info "it seems that all FASTQ files are single-end"
        endedness = "single"
        fqs_dict = Dict(id => Dict("R1" => Dict{String,Vector{String}}()) for id in unique(getindex.(full_fq_ms, 1)))
        for full_fq_m in full_fq_ms
            if !haskey(fqs_dict[full_fq_m[1]]["R1"], full_fq_m[2])
                fqs_dict[full_fq_m[1]]["R1"][full_fq_m[2]] = [full_fq_m[6]]
            else
                push!(fqs_dict[full_fq_m[1]]["R1"][full_fq_m[2]], full_fq_m[6])
            end
        end
    else
        @error "it seems that some FASTQ files are paired-end but others are single-end"
    end

    return [endedness, fqs_dict]
end