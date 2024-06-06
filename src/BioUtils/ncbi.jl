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
	valid_library_types = ("10x", "sf")

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

		if library_type == "sf"
			for acc in accessions
				cmd = `fasterq-dump $acc --threads $nthreads --split-files --outdir $out_dir`
				@info string("running ", cmd, " ...")
				run(cmd)
				@info string("running ", cmd, " done!")
			end
		elseif library_type == "10x"
			for acc in accessions
				cmd = `fasterq-dump $acc --threads $nthreads --split-files --include-technical --outdir $out_dir`
				@info string("running ", cmd, " ...")
				run(cmd)
				@info string("running ", cmd, " done!")
			end
		end
	else
		if library_type == "sf"
			cmd = `fasterq-dump $accession --threads $nthreads --split-files --outdir $out_dir`
			@info string("running ", cmd, " ...")
			run(cmd)
			@info string("running ", cmd, " done!")
		elseif library_type == "10x"
			cmd = `fasterq-dump $accession --threads $nthreads --split-files --include-technical --outdir $out_dir`
			@info string("running ", cmd, " ...")
			run(cmd)
			@info string("running ", cmd, " done!")
		end
	end
end
