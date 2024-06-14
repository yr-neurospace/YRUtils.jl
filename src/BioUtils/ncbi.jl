"""
	ncbi_sra_prefetch(accession::AbstractString, out_dir::AbstractString) -> nothing

A simple wrapper of NCBI SRA `prefetch` command.

`accession` can be either an SRA accession or a file containing many SRA accessions (one per line).
"""
function ncbi_sra_prefetch(accession::AbstractString, out_dir::AbstractString)
	if isempty(accession)
		@error "accession is empty"
	end
	if isempty(out_dir)
		@error "out_dir is empty"
	else
		out_dir = abspath(expanduser(out_dir))
	end

	if isfile(accession)
		accession = abspath(expanduser(accession))
		accessions = open(accession, "r") do io
			x = readlines(io)
			x[.!isempty.(x)]
		end

		for acc in accessions
			cmd = `prefetch $acc -O $out_dir --max-size u`
			@info string("running ", cmd, " ...")
			run(cmd)
			@info string("running ", cmd, " done!")
		end
	else
		cmd = `prefetch $accession -O $out_dir --max-size u`
		@info string("running ", cmd, " ...")
		run(cmd)
		@info string("running ", cmd, " done!")
	end
end

"""
	ncbi_sra_dump(accession::AbstractString, out_dir::AbstractString, library_type::AbstractString, nthreads::Int) -> nothing

A simple wrapper of NCBI SRA `fasterq-dump` command.

`accession` can be either an SRA accession or a file containing many SRA accessions (one per line).

At present, `library_type` only supports `sf` (`--split-files`) and `10x` (`--split-files --include-technical`).
"""
function ncbi_sra_dump(accession::AbstractString, out_dir::AbstractString, library_type::AbstractString, nthreads::Int)
	valid_library_types = ("sf", "10x")

	function gen_sf_cmd(acc::AbstractString, nthreads::Int, out_dir::AbstractString)
		cmd = `fasterq-dump $acc --threads $nthreads --split-files --outdir $out_dir`
		@info string("running ", cmd, " ...")
		run(cmd)
		@info string("running ", cmd, " done!")
	end

	function gen_10x_cmd(acc::AbstractString, nthreads::Int, out_dir::AbstractString)
		cmd = `fasterq-dump $acc --threads $nthreads --split-files --include-technical --outdir $out_dir`
		@info string("running ", cmd, " ...")
		run(cmd)
		@info string("running ", cmd, " done!")
	end

	func_dict = Dict(zip(valid_library_types, (gen_sf_cmd, gen_10x_cmd)))

	if isempty(accession)
		@error "accession is empty"
	end
	if isempty(out_dir)
		@error "out_dir is empty"
	else
		out_dir = abspath(expanduser(out_dir))
	end
	if library_type ∉ valid_library_types
		@error string("supported library types: ", valid_library_types)
	end
	if nthreads < 1
		@error "nthreads must be geater than 0"
	end

	if isfile(accession)
		accession = abspath(expanduser(accession))
		accessions = open(accession, "r") do io
			x = readlines(io)
			x[.!isempty.(x)]
		end

		for acc in accessions
			func_dict[library_type](acc, nthreads, out_dir)
			recur_pigz(out_dir, r"\.(fastq|fq)$"; recursive = true, pigz_options = "", num_jobs = 1)
		end
	else
		func_dict[library_type](accession, nthreads, out_dir)
		recur_pigz(out_dir, r"\.(fastq|fq)$"; recursive = true, pigz_options = "", num_jobs = 1)
	end
end
