"""
    ncbi_sra_prefetch(accession::AbstractString, out_dir::AbstractString)::Nothing

A simple wrapper of NCBI SRA `prefetch` command.

`accession` can be either an SRA accession or a file containing many SRA accessions (one per line).
"""
function ncbi_sra_prefetch(accession::AbstractString, out_dir::AbstractString)::Nothing
    accession = strip(accession)
    out_dir = strip(out_dir)

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
            strip.(x[.!isempty.(x)])
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
    ncbi_sra_dump(accession::AbstractString, out_dir::AbstractString=pwd(), extra_args::AbstractString="", nthreads::Int=6)::Nothing

A simple wrapper of NCBI SRA `fasterq-dump` command.

`accession` can be either an SRA accession or a file containing many SRA accessions (one per line).

When using `fasterq-dump` to dump an SRA file into a number of FASTQ files, 
one thing must be understood - related biological/technical reads stored in different FASTQ files 
are stored in a single spot in an SRA file. Therefore, dumping an SRA file into a number of 
FASTQ files means that splitting/extracting biological/technical reads from spots and putting each of them 
into a correct FASTQ file. So before running `fasterq-dump`, you need to know that which types of reads are 
stored in a given SRA file, how many FASTQ files you want to dump those reads into, and whether you want to 
include technical reads.

You can refer to `https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump` for more detailed explanation
of each parameter. In brief, `fasterq-dump` can work in the following modes:

- `--split-3`: split spots into three types of **biological** reads (two mated, one unmated). For spots having 2 reads, 
    the reads are written into the `*_1.fastq` and `*_2.fastq` files. Unmated reads are placed in `*.fastq` file. 
    This is the default mode. Therefore, by default, `fasterq-dump` can deal with single-end and paired-end reads 
    automatically.

- `--split-spot`: split spots into **biological/technical** reads, and write all of them into a single FASTQ file in order. 
    This mode allows for the output to be redirected to stdout via `--stdout` or `-Z`.

- `--split-files`: split spots into **biological/technical** reads and write each n-th read into a different file. 
        Compared to `--split-3`, the number of files generated depends on the number of read types, which may 
        be more than 3.

- `--concatenate-reads`: spots **including technical reads or not** are NOT split and are written into one FASTQ file as is. 
        This mode allows for the output to be redirected to stdout via `--stdout` or `-Z`.

- `--fasta-unsorted`: this mode is indentical to the `--split-spot` mode, with the only difference being that 
    the original order of the spots and reads is not preserved and it being exlusivly for FASTA. 
    The reason for the existence of this mode is the fact that this mode is faster then 
    the `--split-spot` mode, and does not use temporary files.

As described above, except for the `--split-3` mode, all other modes support whether to keep technical reads or not. 
By default, `--skip-technical` is set, you can use `--include-technical` to enable this.

e.g. 

- For single-end and paired-end reads, just use the default mode `--split-3 --skip-technical`.

- For sequencing technologies, such as 10X transcriptome, technical reads are essential for 
    downstream analysis, and there may be more than 3 FASTQ files. In such a case, use 
    `--split-files --include-technical`.

In addition, parameters `--outdir` (specifying output directory), `--threads` (specifying the number of threads), 
`--seq-defline` (specifying the header contents of sequences, starting with `@` for FASTQ) and `--qual-defline` 
(specifying the header contents of quality sequences, starting with `+` for FASTQ). For available fields, see `fasterq-dump --help`. 
Be careful to choose the correct first character `@`/`+`/`>` based on the desired output (FASTQ/FASTA).
"""
function ncbi_sra_dump(accession::AbstractString, out_dir::AbstractString=pwd(), extra_args::AbstractString="", nthreads::Int=6)::Nothing
    accession = strip(accession)
    out_dir = strip(out_dir)
    extra_args = strip(extra_args)

    if isempty(accession)
        @error "accession is empty"
    end
    if isempty(out_dir)
        @error "out_dir is empty"
    else
        out_dir = abspath(expanduser(out_dir))
    end
    if nthreads < 1
        nthreads = 6
        @warn "nthreads cannot be less than 1, and has been reset to 6"
    end

    if isfile(accession)
        accession = abspath(expanduser(accession))
        accessions = open(accession, "r") do io
            x = readlines(io)
            strip.(x[.!isempty.(x)])
        end

        for acc in accessions
            cmd = if isempty(extra_args)
                `fasterq-dump $acc --threads $nthreads --outdir $out_dir`
            else
                `fasterq-dump $acc --threads $nthreads --outdir $out_dir $extra_args`
            end
            @info string("running ", cmd, " ...")
            run(cmd)
            @info string("running ", cmd, " done!")
            recur_pigz(out_dir, r"\.(fastq|fq)$"; recursive=true, pigz_options="", num_jobs=1)
        end
    else
        acc = accession
        cmd = if isempty(extra_args)
            `fasterq-dump $acc --threads $nthreads --outdir $out_dir`
        else
            `fasterq-dump $acc --threads $nthreads --outdir $out_dir $extra_args`
        end
        @info string("running ", cmd, " ...")
        run(cmd)
        @info string("running ", cmd, " done!")
        recur_pigz(out_dir, r"\.(fastq|fq)$"; recursive=true, pigz_options="", num_jobs=1)
    end
end
