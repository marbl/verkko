# Verkko

Verkko is a hybrid genome assembly pipeline developed for
telomere-to-telomere assembly of PacBio HiFi or Oxford Nanopore Duplex and Oxford Nanopore simplex reads. 
Verkko is Finnish for net, mesh and graph.

Verkko uses Canu to correct remaining errors in the HiFi/duplex reads, builds a
multiplex de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the
graph using GraphAligner, progressively resolves loops and tangles first with
the HiFi reads then with the aligned Oxford Nanopore reads, and finally
creates contig consensus sequences using Canu's consensus module.

## Install:

* Compilation from source requires [GCC 9](https://gcc.gnu.org/) or newer and [Rust 1.66.1](https://www.rust-lang.org/) or newer.  (Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself.)

* Running verkko requires:
  * [Python](https://www.python.org) (v3.5+) 
  * [Snakemake](https://snakemake.readthedocs.io/en/stable/) (>= v7.0, < 8.0.1)
  * [GraphAligner](https://github.com/maickrau/GraphAligner)
  * [MashMap](https://github.com/marbl/MashMap)
  * [Winnowmap](https://github.com/marbl/Winnowmap)
* Running verkko with hi-c data also requires 
  * [Samtools](http://www.htslib.org/)
  * [BWA](https://bio-bwa.sourceforge.net/)
  * [seqtk](https://github.com/lh3/seqtk)
  * [networkx](https://networkx.org/documentation/stable/install.html) python library
    
Installing with a 'package manager' is encouraged:
  * `conda install -c conda-forge -c bioconda -c defaults verkko`
   
or
  * `conda create -n verkko -c conda-forge -c bioconda -c defaults verkko`
  
if you prefer to install verkko in a separate environment. Alternatively, you can download the source for a recent [release](https://github.com/marbl/verkko/releases).

To install Verkko from github (for developers only) run:

    git clone https://github.com/marbl/verkko.git
    cd verkko/src
    git checkout <desired branch> (optional if you want to use a branch for development/compilation and not master)
    make -j32

This will create the folder `verkko/bin` and `verkko/lib/verkko`. You can move the contents of these folders to a central installation location or you can add `verkko/bin` to your path. If any of the dependencies (e.g. GraphAligner, winnowmap, mashmap, etc) are not available in your path you may also symlink them under `verkko/lib/verkko/bin/`.

## Run:

Verkko is implemented as a Snakemake workflow, launched by a wrapper script to parse options
and create a config.yml file.

    verkko -d <work-directory> --hifi <hifi-read-files> [--nano <ont-read-files>]

You can provide any combination of PacBio HiFi/Oxford Nanopore duplex/both to the --hifi parameter. By default, verkko will run the snakemake workflow and all compute on the local machine. Support for SGE, Slurm and LSF (untested) can be enabled with options `--sge`, `--slurm` and `--lsf`, respectively. This will run the snakemake workflow on the local machine but submit all compute to the grid. To launch the both the snakemake workflow and compute on the grid, wrap the verkko command in a shell script and submit using your scheduler.  You may need to set the environment variable VERKKO to the installation directory of Verkko if there are errors that component scripts are not found.

Verkko supports extended phasing using using [rukki](https://github.com/marbl/rukki) using either trio or Hi-C information.

To run in trio mode, you must first generate [merqury](https://github.com/marbl/merqury) hapmer databases and pass them to verkko. Please use git clone to pull the latest versions merqury (see the merqury documentation for details) and make sure that `/path/to/verkko/lib/verkko/bin` is in your path. Then, if you have a SLURM cluster you can run:

    $MERQURY/_submit_build.sh -c 30 maternal.fofn maternal_compress
    $MERQURY/_submit_build.sh -c 30 paternal.fofn paternal_compress
    $MERQURY/_submit_build.sh -c 30 child.fofn    child_compress
_submit_build.sh scripts should not be run from the same directory at the same time since Merqury uses temporary files with common names.
    
if not, you can run

    meryl count compress k=30 threads=XX memory=YY maternal.*fastq.gz output maternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY paternal.*fastq.gz output paternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY    child.*fastq.gz output    child_compress.k30.meryl

replacing XX and YY with the threads and memory you want meryl to use. Once you have the databases, run:

    $MERQURY/trio/hapmers.sh \
      maternal_compress.k30.meryl \
      paternal_compress.k30.meryl \
         child_compress.k30.meryl

    verkko -d asm \
      --hifi hifi/*.fastq.gz \
      --nano  ont/*.fastq.gz \
      --hap-kmers maternal_compress.k30.hapmer.meryl \
                  paternal_compress.k30.hapmer.meryl \
                  trio

Make sure to count k-mers in compressed space. Child data is optional, in this case use `maternal_compress.k30.only.meryl` and  `paternal_compress.k30.only.meryl` in the verkko command above.

To run in Hi-C mode, reads should be provided using the --hic1 and --hic2 options. For example:

    verkko -d asm \
      --hifi hifi/*.fastq.gz \
      --nano ont/*.fastq.gz \
      --hic1 hic/*R1*fastq.gz  \
      --hic2 hic/*R2*fastq.gz

Hi-C integration was tested mostly on human and primate genomes. Please see the --rdna_scaff_ref, --rdna-tangle, --uneven-depth and --haplo-divergence options if you want to assemble something distant from human and/or have uneven coverage. If you encounter issues or have questions about appropriate parameters, please open an [issue](https://github.com/marbl/verkko/issues).

For slurm/sge cluster runs verkko uses different cpu/memory/time options for different parts of the pipeline. Usually user is not supposed to change them, however advanced tuning is possible with --<stage_code>-run options.
<details>
<summary>Here we list those options and briefly describe of corresponding verkko pipeline stages and affected scripts </summary>
 
    --sto-run                Creating read storage, 0-correction/buildStore.sh
    --ovb-run                HiFi/Duplex read overlapping,  0-correction/matchchains-index.sh & 0-correction/overlap-jobs/*.sh
    --ovs-run                Combining overlap info, 0-correction/combineOverlapsConfigure.sh & 0-correction/combineOverlaps.sh
    --red-run                Read correction, 0-correction/configureFindErrors.sh, 0-correction/find-errors-jobs/*.sh & 0-correction/fixErrors.sh
    --mbg-run                Multiplex de Bruijn graph construction, 1-buildGraph/buildGraph.sh
    --utg-run                Initial graph simplification, 2-processGraph/processGraph.sh
    --spl-run                ONT reads splitting, 3-alignTips/splitONT.sh & 3-align/splitONT.sh
    --ali-run                ONT reads alignment, 3-align/aligned*.sh & 3-alignTips/aligned*.sh
    --pop-run                ONT-based graph simplification, 4-processONT/processONT.sh
    --utp-run                Final graph simplification(tip clipping), 5-untip/untip.sh
    --lay-run                Contigs preprocessing for consensus, 6-layoutContigs/createLayoutInputs.sh & 6-layoutContigs/createLayout.sh
    --sub-run                Extraction of the read subset for consensus,  7-consensus/extractONT.sh
    --par-run                Reads preprocessing for consensus, 7-consensus/buildPackages.sh
    --cns-run                Reads consenus, 7-consensus/packages/part*.sh
    --ahc-run                HiC alignment, 8-hicPipeline/align_bwa*.sh
    --fhc-run                All scripts in HiC pipeline other than alignment,  8-hicPipeline/*.sh

Default values can be found in verkko bash script, i.e.  `grep par_ bin/verkko`. Values used for verkko runs are listed in verkko.yml in the run directory.
</details> 

You can pass through snakemake options to restrict CPU/memory/cluster resources by adding the `--snakeopts` option to verkko. For example, `--snakeopts "--dry-run"` will print what jobs will run while `--snakeopts "--cores 1000"` would restrict grid runs to at most 1000 cores across all submited jobs.

To test your installation we have an E. coli K12 dataset available. 

    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o hifi.fastq.gz
    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o ont.fastq.gz
    verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz

The final assembly result is under `asm/assembly.fasta`. The final graph (in homopolymer-compressed space) is under `asm/assembly.homopolymer-compressed.gfa` along with coverage files in `asm/assembly*csv`. If you provided phasing information, you will also have `asm/assembly.haplotype[12].fasta`. You can find intermediate graphs and coverage files under `asm/*/unitig-*gfa` and `asm/*/unitig-*csv`.

## Citations:
 - Rautiainen M, Nurk S, Walenz BP, Logsdon GA, Porubsky D, Rhie A, Eichler EE, Phillippy AM, Koren S. [Telomere-to-telomere assembly of diploid chromosomes with Verkko](https://doi.org/10.1038/s41587-023-01662-6). Nat Biotech. (2023). `doi:10.1038/s41587-023-01662-6`
