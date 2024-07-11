"""
	rename_10x_fastq_name(
	dir::AbstractString,
	sra_run_table::AbstractString,
	sample_name_func::Function;
	rna_read_scheme::Dict{Int64, String} = Dict(8 => "I1", 28 => "R1", 90 => "R2"),
	atac_read_scheme::Dict{Int64, String} = Dict(8 => "I1", 16 => "I2", 50 => "R1", 49 => "R2"),
	skip::Bool = false,
	out_dir::Vector{String} = ["rna", "atac"],
	pattern::Regex = r"") -> Tuple{String, String}

Rename 10X scRNA-seq/scATAC-seq FASTQ files dumped from NCBI GEO SRA files to make them suitable to run the 10X `cellranger` pipeline.

`sra_run_table` is a metadata file downloaded from NCBI GEO.

`sample_name_func`, accepting a single `DataFrameRow` argument, is used to get the sample name from each row of `sra_run_table`.

`out_dir` specifies the output directories of renamed FASTQ files for scRNA-seq and scATAC-seq respectively.

The first is for scRNA-seq, and the second is for scATAC-seq.

If `skip` is `true`, then items having no corresponding FASTQ files will be skipped.

Finally, files named "rename.sh" and "libraries.csv" will be created in the directory `dir`, and then you can check whether they are reasonable.

If OK, just run "rename.sh" in your Linux terminal.

# Examples
```julia-repl
julia> function sample_name_func(x::DataFrameRow)
		   x.donor_id
	   end
```
"""
function rename_10x_fastq_name(
	dir::AbstractString,
	sra_run_table::AbstractString,
	sample_name_func::Function;
	rna_read_scheme::Dict{Int64, String} = Dict(8 => "I1", 28 => "R1", 90 => "R2"),
	atac_read_scheme::Dict{Int64, String} = Dict(8 => "I1", 16 => "I2", 50 => "R1", 49 => "R2"),
	skip::Bool = false,
	out_dir::Vector{String} = ["rna", "atac"],
	pattern::Regex = r"")

	if isempty(dir)
		@error "dir is empty"
	else
		dir = abspath(expanduser(dir))
	end

	if !isfile(sra_run_table)
		@error "sra_run_table is invalid"
	else
		sra_run_table = abspath(expanduser(sra_run_table))
	end

	if isempty(out_dir) || length(out_dir) != 2 || any(isempty.(out_dir))
		@error "out_dir is invalid"
	else
		out_dir = joinpath.(dir, out_dir)
	end

	if pattern == r""
		pattern = r"\.(fastq|fq|fastq\.gz|fq\.gz)$"
	end

	cmds = string.("mkdir -p ", out_dir)
	lib_csv = ["fastqs,sample,library_type"]

	sra_run_table_df = CSV.read(sra_run_table, DataFrame; delim = ',')
	if isempty(sra_run_table_df)
		@error "the content read from sra_run_table is empty"
	end

	all_fq_files = list_files(dir, pattern; recursive = false, full_name = true)
	if isempty(all_fq_files)
		@error string("no FASTQ files found in the directory: ", dir)
	end

	for r in eachrow(sra_run_table_df)
		cmds = [cmds; "\n"]
		r_fq_files = all_fq_files[.!isnothing.(match.(Regex(string("^", r.Run)), basename.(all_fq_files)))]
		if isempty(r_fq_files)
			if skip
				continue
			else
				@error string(r.Run, " didn't match any FASTQ files")
			end
		end

		for r_fq_file in r_fq_files
			lens = FASTQReader(GzipDecompressorStream(open(r_fq_file))) do reader
				count = 0
				lens = Array{Int, 1}(undef, 10000)
				for record in reader
					count += 1
					if count > 10000
						break
					end
					lens[count] = length(sequence(record))
				end
				return unique(lens)
			end

			if length(lens) != 1
				@error string("the lengthes of the first 1e4 reads in ", r_fq_file, " are not consistent")
			end

			if r."Assay Type" == "RNA-Seq"
				if haskey(rna_read_scheme, lens[1])
					cmds = [
						cmds;
						string("mv ", r_fq_file, " ",
							joinpath(out_dir[1], string(sample_name_func(r), "_S1_L001_", rna_read_scheme[lens[1]], "_001.fastq.gz")))
					]
					lib_csv = [
						lib_csv;
						string(out_dir[1], ",", sample_name_func(r), ",", "Gene Expression")
					]
				else
					@error string("no matched key for ", lens[1], " in RNA-seq scheme")
				end
			elseif r."Assay Type" == "ATAC-seq"
				if haskey(atac_read_scheme, lens[1])
					cmds = [
						cmds;
						string("mv ", r_fq_file, " ",
							joinpath(out_dir[2], string(sample_name_func(r), "_S1_L001_", atac_read_scheme[lens[1]], "_001.fastq.gz")))
					]
					lib_csv = [
						lib_csv;
						string(out_dir[2], ",", sample_name_func(r), ",", "Chromatin Accessibility")
					]
				else
					@error string("no matched key for ", lens[1], " in ATAC-seq scheme")
				end
			else
				@error string("unsupported assay type: ", r."Assay Type")
			end
		end
	end

	rename_sh_file = joinpath(dir, "rename.sh")
	open(rename_sh_file, "w") do io
		println.(io, cmds)
	end

	lib_csv_file = joinpath(dir, "libraries.csv")
	open(lib_csv_file, "w") do io
		println.(io, unique(lib_csv))
	end

	return (rename_sh_file, lib_csv_file)
end

"""
	run_10x_cellranger_arc_count(
	library::AbstractString,
	reference::AbstractString,
	localcores::Int = 0,
	localmem::Int = 0,
	working_dir::AbstractString = "",
	log_file::AbstractString = "run_10x_cellranger_arc_count.log",
	cluster_header::AbstractString = "#PBS -N 10x_cellranger_arc_count\n#PBS -l nodes=1:ppn=60\n#PBS -l walltime=7200:00:00\n#PBS -q mem3T\n#PBS -V") -> String

Generate the shell script file for running 10X `cellranger-arc count` in batch.
"""
function run_10x_cellranger_arc_count(
	library::AbstractString,
	reference::AbstractString,
	localcores::Int = 0,
	localmem::Int = 0,
	working_dir::AbstractString = "",
	log_file::AbstractString = "run_10x_cellranger_arc_count.log",
	cluster_header::AbstractString = "#PBS -N 10x_cellranger_arc_count\n#PBS -l nodes=1:ppn=24\n#PBS -l walltime=500:00:00\n#PBS -V")

	if isempty(library)
		@error "library is empty"
	end
	if isempty(reference)
		@error "reference is empty"
	end
	if localcores < 1
		@info "won't set the argument --localcores"
	end
	if localmem < 1
		@info "won't set the argument --localmem"
	end

	output_sh_file = joinpath(dirname(library), "run_10x_cellranger_arc_count.sh")
	lib_df = CSV.read(library, DataFrame; delim = ',')
	samples = unique(lib_df[:, "sample"])
	open(output_sh_file, "w") do io
		if !isempty(cluster_header)
			println(io, cluster_header)
		end
		if !isempty(working_dir)
			println(io, "\n", "cd ", working_dir)
		end

		for sample in samples
			cmd = string(
				"cellranger-arc count --id=$sample --reference=$reference --libraries=$library",
				if localcores > 0
					" --localcores=$localcores"
				else
					""
				end,
				if localmem > 0
					" --localmem=$localmem"
				else
					""
				end,
				if !isempty(log_file)
					" &>> $log_file"
				else
					""
				end)
			println(io, "\n", cmd)
		end
	end

	return output_sh_file
end

"""
	run_10x_cellranger_count(
	library::AbstractString,
	reference::AbstractString,
	localcores::Int = 0,
	localmem::Int = 0,
	create_bam::Bool = false,
	working_dir::AbstractString = "",
	log_file::AbstractString = "run_10x_cellranger_count.log",
	cluster_header::AbstractString = "#PBS -N 10x_cellranger_count\n#PBS -l nodes=1:ppn=24\n#PBS -l walltime=500:00:00\n#PBS -V") -> String

Generate the shell script file for running 10X `cellranger count` in batch.
"""
function run_10x_cellranger_count(
	library::AbstractString,
	reference::AbstractString,
	localcores::Int = 0,
	localmem::Int = 0,
	create_bam::Bool = false,
	working_dir::AbstractString = "",
	log_file::AbstractString = "run_10x_cellranger_count.log",
	cluster_header::AbstractString = "#PBS -N 10x_cellranger_count\n#PBS -l nodes=1:ppn=24\n#PBS -l walltime=500:00:00\n#PBS -V")

	if isempty(library)
		@error "library is empty"
	end
	if isempty(reference)
		@error "reference is empty"
	end
	if localcores < 1
		@info "won't set the argument --localcores"
	end
	if localmem < 1
		@info "won't set the argument --localmem"
	end

	output_sh_file = joinpath(dirname(library), "run_10x_cellranger_arc_count.sh")
	lib_df = CSV.read(library, DataFrame; delim = ',')
	open(output_sh_file, "w") do io
		if !isempty(cluster_header)
			println(io, cluster_header)
		end
		if !isempty(working_dir)
			println(io, "\n", "cd ", working_dir)
		end

		for r in eachrow(lib_df)
			cmd = string(
				"cellranger count --id=$(r.sample) --transcriptome=$reference --fastqs=$(r.fastqs) --sample=$(r.sample)",
				if create_bam
					" --create-bam=true"
				else
					""
				end,
				if localcores > 0
					" --localcores=$localcores"
				else
					""
				end,
				if localmem > 0
					" --localmem=$localmem"
				else
					""
				end,
				if !isempty(log_file)
					" &>> $log_file"
				else
					""
				end)
			println(io, "\n", cmd)
		end
	end

	return output_sh_file
end