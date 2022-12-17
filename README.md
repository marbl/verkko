# Verkko

Verkko is a hybrid genome assembly pipeline developed for
telomere-to-telomere assembly of PacBio HiFi and Oxford Nanopore reads.
Verkko is Finnish for net, mesh and graph.

Verkko uses Canu to correct remaining errors in the HiFi reads, builds a
multiplex de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the
graph using GraphAligner, progressively resolves loops and tangles first with
the HiFi reads then with the aligned Oxford Nanopore reads, and finally
creates contig consensus sequences using Canu's consensus module.

**Note: Verkko is a work in progress and currently available as a beta
release. Expect to encounter instability, but please feel free to submit
an issue if you encounter any problems.**

## Install:

* Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself. Compilation from source requires GCC 7 or newer and Rust 1.58 or newer.

* Running verkko requires both [MBG](https://github.com/maickrau/MBG) and [GraphAligner](https://github.com/maickrau/GraphAligner) to be installed. Verkko also requires [Snakemake](https://snakemake.readthedocs.io/en/stable/) (v7.0+) and a recent [Python](https://www.python.org) (v3.5+). 

Installing with a 'package manager' is encouraged:
  * `conda install -c conda-forge -c bioconda -c defaults verkko`
   
or
  * `conda create -n verkko -c conda-forge -c bioconda -c defaults verkko`
  
if you prefer to install verkko in a separate environment. Alternatively, you can download the source for a recent [release](https://github.com/maickrau/verkko/releases).

To install Verkko from github (for developers only) run:

    git clone https://github.com/maickrau/verkko.git
    cd verkko/src
    git submodule init && git submodule update
    make -j32

This will create the folder `verkko/bin` and `verkko/lib/verkko`. You can move the contents of these folders to a central installation location or you can add `verkko/bin` to your path. If GraphAligner or MBG are not available in your path you may also symlink them under `verkko/lib/verkko/bin/`. Make sure you are using the latest tip of MBG/GraphAligner not a conda install in this case.

## Run:

Verkko is implemented as a Snakemake workflow, launched by a wrapper script to parse options
and create a config.yml file.

    verkko -d <work-directory> --hifi <hifi-read-files> [--nano <ont-read-files>]

By default, verkko will run the snakemake workflow and all compute on the local machine. Support for SGE, Slurm and LSF (untested) can be enabled with options `--sge`, `--slurm` and `--lsf`, respectively. This will run the snakemake workflow on the local machine but submit all compute to the grid. To launch the both the snakemake workflow and compute on the grid, wrap the verkko command in a shell script and submit using your scheduler.  You may need to set the environment variable VERKKO to the installation directory of Verkko if there are errors that component scripts are not found.

Verkko supports trio-based phasing using using [rukki](https://github.com/marbl/rukki). To run in this mode, you must first generate [merqury](https://github.com/marbl/merqury) hapmer databases and pass them to verkko. Please use git clone to pull the latest versions merqury (see the merqury documentation for details) and make sure that `/path/to/verkko/lib/verkko/bin` is in your path. Then, if you have a SLURM cluster you can run:

    $MERQURY/_submit_build.sh -c 30 maternal.fofn maternal_compress
    $MERQURY/_submit_build.sh -c 30 paternal.fofn paternal_compress
    $MERQURY/_submit_build.sh -c 30 child.fofn    child_compress
    
if not, you can run

    meryl count compress k=30 threads=XX memory=YY maternal.*fastq.gz output maternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY paternal.*fastq.gz output paternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY child.*fastq.gz output paternal_compress.k30.meryl

replacing XX and YY with the threads and memory you want meryl to use. Once you have the databases, run:

    $MERQURY/trio/hapmers.sh maternal_compress.k30.meryl paternal_compress.k30.meryl child_compress.k30.meryl
    verkko -d asm --hifi hifi/*.fastq.gz --nano ont/*.fastq.gz --hap-kmers maternal_compress.k30.hapmer.meryl paternal_compress.k30.hapmer.meryl trio

Make sure to count k-mers in compressed space. Child data is optional, in this case use `maternal_compress.k30.only.meryl` and  `paternal_compress.k30.only.meryl` in the verkko command above. Preliminary support is available for read sets binned by haplotype from another method, such as [PGAS](https://github.com/daewoooo/SaaRclust) and Strand-Seq or [DipAsm](https://github.com/shilpagarg/DipAsm) and Hi-C. In these cases, make sure the phase blocks are chromosome-scale and consistent within each chromosome. You can build merqury DBs as above and specify them along with either `hic` or `strandseq` instead of `trio` to verkko instead.

You can pass through snakemake options to restrict CPU/memory/cluster resources by adding the `--snakeopts` option to verkko. For example, `--snakeopts "--dry-run"` will print what jobs will run while `--snakeopts "--cores 1000"` would restrict grid runs to at most 1000 cores across all submited jobs.

To test your installation we have an E. coli K12 dataset available. 

    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o hifi.fastq.gz
    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o ont.fastq.gz
    verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz

The final assembly result is under `asm/assembly.fasta`. The final graph (in homopolymer-compressed space) is under `asm/assembly.homopolymer-compressed.gfa` along with coverage files in `asm/assembly*csv`. If you provided phasing information, you will also have `asm/assembly.haplotype[12].fasta`. You can find intermediate graphs and coverage files under `asm/*/unitig-*gfa` and `asm/*/unitig-*csv`.

## Citations:
 - Rautiainen M, Nurk S, Walenz BP, Logsdon GA, Porubsky D, Rhie A, Eichler EE, Phillippy AM, Koren S. [Verkko: telomere-to-telomere assembly of diploid chromosomes](https://doi.org/10.1101/2022.06.24.497523). bioRxiv. (2022). `doi:10.1101/2022.06.24.497523`
