"""
    fastqc(files::Vector{String}, outdir::AbstractString=pwd();
    fastqc_path::AbstractString="", fastqc_options::AbstractString="--threads 1",
    multiqc_path::AbstractString="", multiqc_options::AbstractString="--zip-data-dir",
    kwargs...)::Vector{String}

Call `fastqc` in Julia on FASTQ `files` followed by calling `multiqc` on the results of `fastqc`.

The path to `fastqc`/`multiqc` will be found using `which fastqc`/`which multiqc` if `fastqc_path = ""`/`multiqc_path = ""`. 

Other `fastqc`/`multiqc` options can be given by `fastqc_options`/`multiqc_options`.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function fastqc(files::Vector{String}, outdir::AbstractString=pwd();
    fastqc_path::AbstractString="", fastqc_options::AbstractString="--threads 1",
    multiqc_path::AbstractString="", multiqc_options::AbstractString="--zip-data-dir",
    kwargs...)::Vector{String}
    if isempty(files) || isempty(strip(outdir))
        @error "both files and outdir cannot be empty"
    else
        files = abspath.(expanduser.(files))
        outdir = abspath(expanduser(outdir))
    end

    if isempty(strip(fastqc_path))
        fastqc_path = find_cmd("fastqc"; return_nothing=false)
    end
    cmd_valid(Cmd(string.([fastqc_path, "--version"])); return_false=false)

    if isempty(strip(multiqc_path))
        multiqc_path = find_cmd("multiqc"; return_nothing=false)
    end
    cmd_valid(Cmd(string.([multiqc_path, "--version"])); return_false=false)

    fastqc_cmd_vec = [fastqc_path, "--outdir", outdir]
    fastqc_options = strip(fastqc_options)
    if !isempty(fastqc_options)
        fastqc_cmd_vec = [fastqc_cmd_vec; Base.shell_split(fastqc_options)]
    end

    @info "start running fastqc ..."
    para_cmds(files; kwargs...) do x
        stdout_io_bf = IOBuffer()
        stderr_io_bf = IOBuffer()
        cmd = Cmd(string.([fastqc_cmd_vec; x]))
        run(pipeline(cmd; stdout=stdout_io_bf, stderr=stderr_io_bf); wait=true)
        stdout_stderr_str = string(String(take!(stdout_io_bf)),
            String(take!(stderr_io_bf)), "\n\n")
        print(stdout_stderr_str)
    end

    multiqc_cmd_vec = [multiqc_path, "--outdir", outdir, outdir]
    multiqc_options = strip(multiqc_options)
    if !isempty(multiqc_options)
        multiqc_cmd_vec = [multiqc_cmd_vec; Base.shell_split(multiqc_options)]
    end

    @info "start running multiqc ..."
    run(Cmd(string.(multiqc_cmd_vec)); wait=true)

    files
end

"""
    trimgalore(dict::Union{Dict{String,Dict{String,Dict{String,Vector{String}}}},Dict{String,Dict{String,Vector{String}}}}, read_type::AbstractString="paired", outdir::AbstractString=pwd();
    trimgalore_path::AbstractString="", trimgalore_options::AbstractString="--cores 4 --phred33 --quality 20 --length 30 --trim-n", kwargs...)::Vector

Call `trim_galore` in Julia over FASTQ files given in `dict`.

`read_type` can only be either `paired` or `single`.

`outdir` must have already been created in advance.

The `dict` for paired end is of type `Dict{String, Dict{String, Dict{String, Vector{String}}}}`.
i.e.

- The first `String` is sample IDs.

- The second `String` is biological replicate IDs.

- The third `String` is read type IDs (must exactly be `R1` and `R2`).

- The fourth `Vector{String}` is a list of FASTQ files which may contain technical replicates.

Note: the order of technical replicates must be exactly same between `R1` and `R2`.

The `dict` for single end is of type `Dict{String, Dict{String, Vector{String}}}`.
In comparison with the paired end, the only difference is that the single end lacks the read type IDs.
i.e. the technical replicates directly follow the biological replicate IDs instead of the read type IDs.

For how to generate `dict`, see `auto_detect_fastq_read_type`.

The path to `trim_galore` will be found using `which trim_galore` if `trimgalore_path = ""`. 

Other `trim_galore` options can be given by `trimgalore_options`.

All other keyword arguments will be passed to `para_cmds` via `kwargs`.
"""
function trimgalore(dict::Union{Dict{String,Dict{String,Dict{String,Vector{String}}}},Dict{String,Dict{String,Vector{String}}}}, read_type::AbstractString="paired", outdir::AbstractString=pwd();
    trimgalore_path::AbstractString="", trimgalore_options::AbstractString="--cores 4 --phred33 --quality 20 --length 30 --trim-n", kwargs...)::Vector
    if isempty(strip(outdir))
        @error "outdir cannot be empty"
    else
        outdir = abspath(expanduser(outdir))
    end

    if !(read_type in ("paired", "single"))
        @error "read_type can only be either paired or single"
    end

    if isempty(strip(trimgalore_path))
        trimgalore_path = find_cmd("trim_galore")
    end

    cmd_valid(Cmd(string.([trimgalore_path, "--version"])); return_false=false)

    cutadapt_path = find_cmd("cutadapt")
    pigz_path = find_cmd("pigz")
    cmd_valid(Cmd(string.([cutadapt_path, "--version"])); return_false=false)
    cmd_valid(Cmd(string.([pigz_path, "--version"])); return_false=false)

    trimgalore_options = strip(trimgalore_options)
    cmd_vec = [trimgalore_path]
    if !isempty(trimgalore_options)
        cmd_vec = [cmd_vec; Base.shell_split(trimgalore_options)]
    end
    cmd_vec = [cmd_vec; "--output_dir"; outdir]

    if read_type == "paired"
        all_r1_r2_files = []
        for (_, rep_dict) in dict
            for (_, read_type_dict) in rep_dict
                all_r1_r2_files = [all_r1_r2_files; collect(zip(read_type_dict["R1"], read_type_dict["R2"]))]
            end
        end

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
    end

    if read_type == "single"
        all_r1_files = []
        for (_, rep_dict) in dict
            for (_, rep_vec) in rep_dict
                all_r1_files = [all_r1_files; rep_vec]
            end
        end

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
    fq_suffix_pattern = r"\.(fastq|fq|fastq\.gz|fq\.gz)$"
    trimmed_fq_files = list_files(outdir, fq_suffix_pattern; recursive=false, full_name=true)
    if isempty(trimmed_fq_files)
        @error "do not find any trimmed FASTQ files"
    end

    renamed_trimmed_fq_files = replace.(trimmed_fq_files,
        Regex(string("(_val_[12]|_trimmed)(?=", fq_suffix_pattern.pattern, ")")) => "")
    for (s, t) in zip(trimmed_fq_files, renamed_trimmed_fq_files)
        @info string("move ", s, " to ", t, " ...")
        mv(s, t; force=false)
    end

    return if read_type == "paired"
        all_r1_r2_files
    else
        all_r1_files
    end
end

"""
    auto_detect_fastq_read_type(files::Vector{String}; allow_mix::Bool=false)::Dict{String,Dict{String,Any}}

Detect FASTQ read types and organize FASTQ files in an ordered way.

Valid FASTQ file formats:

- Single-end: `ID_repN[_partN].(fastq|fq).gz`

- Paired-end: `ID_repN[_partN].R(1|2).(fastq|fq).gz`

`ID` is the sample name, which can only contain `[a-zA-Z0-9]` and does **NOT** start with `[0-9]`.

`repN` means the `N`th biological replicate.

`partN` means the `N`th technical replicate. All technical replicates with the same `ID_repN` should be merged into a single FASTQ file before running downstream analyses.

`N` can only contain `[0-9]`.
"""
function auto_detect_fastq_read_type(files::Vector{String}; allow_mix::Bool=false)::Dict{String,Dict{String,Any}}
    paired_pattern = r"(?<ID>^[a-zA-Z]+[a-zA-Z0-9]*)_(?<BioRep>rep[0-9]+)(_(?<TechRep>part[0-9]+))*\.(?<ReadType>R[12])\.(?<Ext>(fastq|fq)\.gz$)"
    single_pattern = r"(?<ID>^[a-zA-Z]+[a-zA-Z0-9]*)_(?<BioRep>rep[0-9]+)(_(?<TechRep>part[0-9]+))*\.(?<Ext>(fastq|fq)\.gz$)"

    if isempty(files)
        @error "files cannot be empty"
    end

    paired_dict = Dict("file" => String[], "basename" => String[],
        "Sample_Group" => String[],
        "BioRep_Group" => String[],
        "TechRep_Group" => String[],
        "ID" => String[], "BioRep" => String[],
        "TechRep" => Union{Nothing,String}[],
        "ReadType" => String[], "Ext" => String[])
    for file in files
        m = match(paired_pattern, basename(file))
        if !isnothing(m)
            push!(paired_dict["file"], file)
            push!(paired_dict["basename"], basename(file))
            push!(paired_dict["Sample_Group"], m["ID"])
            push!(paired_dict["BioRep_Group"], string(m["ID"], "_", m["BioRep"]))
            push!(paired_dict["TechRep_Group"], string(m["ID"], "_", m["BioRep"], if isnothing(m["TechRep"])
                ""
            else
                string("_", m["TechRep"])
            end))
            push!(paired_dict["ID"], m["ID"])
            push!(paired_dict["BioRep"], m["BioRep"])
            push!(paired_dict["TechRep"], m["TechRep"])
            push!(paired_dict["ReadType"], m["ReadType"])
            push!(paired_dict["Ext"], m["Ext"])
        end
    end
    paired_df = DataFrame(paired_dict)

    single_dict = Dict("file" => String[], "basename" => String[],
        "Sample_Group" => String[],
        "BioRep_Group" => String[],
        "TechRep_Group" => String[],
        "ID" => String[], "BioRep" => String[],
        "TechRep" => Union{Nothing,String}[],
        "Ext" => String[])
    for file in files
        m = match(single_pattern, basename(file))
        if !isnothing(m)
            push!(single_dict["file"], file)
            push!(single_dict["basename"], basename(file))
            push!(single_dict["Sample_Group"], m["ID"])
            push!(single_dict["BioRep_Group"], string(m["ID"], "_", m["BioRep"]))
            push!(single_dict["TechRep_Group"], string(m["ID"], "_", m["BioRep"], if isnothing(m["TechRep"])
                ""
            else
                string("_", m["TechRep"])
            end))
            push!(single_dict["ID"], m["ID"])
            push!(single_dict["BioRep"], m["BioRep"])
            push!(single_dict["TechRep"], m["TechRep"])
            push!(single_dict["Ext"], m["Ext"])
        end
    end
    single_df = DataFrame(single_dict)

    if (nrow(single_df) + nrow(paired_df)) != length(files)
        @error string("the sum of the number of paired-end files and the number of single-end files is not equal to the total number of files.\n",
            "This means that there were some files not parsed (i.e. their naming format is incorrect)")
    end

    if !allow_mix && !isempty(paired_df) && !isempty(single_df)
        @error string("in all the files, ",
            nrow(paired_df), " were parsed as paired-end, and ",
            nrow(single_df), " were parsed as single-end.\n",
            "If your files do contain two types of files, set allow_mix = true instead")
    end

    if !isempty(single_df)
        single_res = Dict(key["ID"] => Dict(r["BioRep"] => String[] for r in eachrow(subdf)) for (key, subdf) in pairs(groupby(unique(select(single_df, :ID, :BioRep)), :ID)))
        for r in eachrow(single_df)
            push!(single_res[r["ID"]][r["BioRep"]], r["file"])
        end
    end

    if !isempty(paired_df)
        split_paired_dict = Dict(key["ReadType"] => subdf for (key, subdf) in pairs(groupby(paired_df, :ReadType)))
        joined_paired_df = outerjoin(split_paired_dict["R1"], split_paired_dict["R2"];
            on=:TechRep_Group, source=:source, validate=(true, true), renamecols="_R1" => "_R2")

        if !all(joined_paired_df[!, :source] .== "both")
            @error "there are files only containing either R1 or R2"
        end

        paired_res = Dict(key["ID"] => Dict(r["BioRep"] => Dict("R1" => String[], "R2" => String[]) for r in eachrow(subdf)) for (key, subdf) in pairs(groupby(unique(select(paired_df, :ID, :BioRep)), :ID)))
        for r in eachrow(joined_paired_df)
            push!(paired_res[r["ID_R1"]][r["BioRep_R1"]][r["ReadType_R1"]], r["file_R1"])
            push!(paired_res[r["ID_R2"]][r["BioRep_R2"]][r["ReadType_R2"]], r["file_R2"])
        end
    end

    Dict("paired" => Dict("data_frame" => paired_df,
            "dict" => if !isempty(paired_df)
                paired_res
            else
                Dict()
            end,
            "status" => if !isempty(paired_df)
                "yes"
            else
                "no"
            end),
        "single" => Dict("data_frame" => single_df,
            "dict" => if !isempty(single_df)
                single_res
            else
                Dict()
            end,
            "status" => if !isempty(single_df)
                "yes"
            else
                "no"
            end))
end