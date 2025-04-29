# Verkko

Verkko is a hybrid genome assembly pipeline developed for
telomere-to-telomere assembly of accurate long reads (PacBio HiFi, Oxford Nanopore Duplex, [HERRO](https://github.com/lbcb-sci/herro) corrected Oxford Nanopore Simplex) and Oxford Nanopore ultra-long reads.
Verkko is Finnish for net, mesh and graph.

Verkko uses Canu to correct remaining errors in the reads, builds a
multiplex de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the
graph using GraphAligner, progressively resolves loops and tangles first with
the HiFi reads then with the aligned Oxford Nanopore reads, and finally
creates contig consensus sequences using Canu's consensus module.

## Table of contents
- [Install](#install)
- [Getting Started](#getting-started)
- [Outputs](#outputs)
- [Test data](#test-data)

## Install:

Installing with a 'package manager' is recommended:
  * `conda create -n verkko -c conda-forge -c bioconda -c defaults verkko`

Alternatively, you can download and compile the source for a recent [release](https://github.com/marbl/verkko/releases).

<details>
<summary><b>Compile from source</b></summary>
 
* Compilation from source requires:
  * [GCC 9](https://gcc.gnu.org/) or newer
  * [Rust 1.74](https://www.rust-lang.org/) or newer.  
 
(Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself.)

* Running verkko requires:
  * [Python](https://www.python.org) (v3.5+)
  * [Snakemake](https://snakemake.readthedocs.io/en/stable/) (>= v7.0, < 8.0.1)
  * [GraphAligner](https://github.com/maickrau/GraphAligner)
  * [MashMap](https://github.com/marbl/MashMap)
  * [Winnowmap](https://github.com/marbl/Winnowmap)
* Running verkko with hi-c/porec data also requires
  * [Samtools](http://www.htslib.org/)
  * [BWA](https://bio-bwa.sourceforge.net/)
  * [Minimap2](https://github.com/lh3/minimap2)
  * [seqtk](https://github.com/lh3/seqtk)
  * [networkx](https://networkx.org/documentation/stable/install.html) python library

To install an unreleased version of Verkko from github (for development) run:

    git clone https://github.com/marbl/verkko.git
    cd verkko/src
    git checkout <desired branch> (optional if you want to use a branch for development/compilation and not master)
    make -j32

This will create the folder `verkko/bin` and `verkko/lib/verkko`. You can move the contents of these folders to a central installation location or you can add `verkko/bin` to your path. If any of the dependencies (e.g. GraphAligner, winnowmap, mashmap, etc) are not available in your path you may also symlink them under `verkko/lib/verkko/bin/`.
</details>

## Getting started:

Verkko is implemented as a Snakemake workflow, launched by a wrapper script to parse options
and create a `verkko.yml` file.

    verkko -d <work-directory> --hifi <hifi-read-files> [--nano <ont-read-files>]

Run `verkko` with no options will list all available options with brief descriptions. At the minimum verkko requires high-accuracy long reads, provided with the `--hifi` option. You can provide any combination of PacBio HiFi/Oxford Nanopore duplex/both to the `--hifi` parameter. However, we strongly recommend including some ultra-long sequence data using the `--nano` parameter and phasing information (see below). For HERRO corrected reads, provide the corrected reads with the `--hifi` option and the uncorrected reads as `--nano`. The output of verkko will be phased scaffolds. Note that no attempt is made to generate a primary or pseudo-haplotype assembly.

### Phasing:
Verkko supports extended phasing using using [rukki](https://github.com/marbl/rukki) using either trio or Hi-C information.

To run in trio mode, you must first generate [merqury](https://github.com/marbl/merqury) hapmer databases and pass them to verkko.
<details>
<summary><b>Build meryl DBs</b></summary>
Please use git clone to pull the latest versions merqury (see the merqury documentation for details). Then, if you have a SLURM cluster you can run:

    # assumes you have maternal/paternal folders
    # each containing a fofn of sequence inputs named [mp]aternal.fofn
    # and a top level folder with a child.fofn specifying F1 sequence inputs
    cd maternal
    $MERQURY/_submit_build.sh -c 30 maternal.fofn maternal_compress
    cd ../paternal
    $MERQURY/_submit_build.sh -c 30 paternal.fofn paternal_compress
    cd ../
    $MERQURY/_submit_build.sh -c 30 child.fofn    child_compress
    ln -s maternal/maternal_compress.k30.meryl
    ln -s paternal/paternal_compress.k30.meryl

without a grid, you can run

    meryl count compress k=30 threads=XX memory=YY maternal.*fastq.gz output maternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY paternal.*fastq.gz output paternal_compress.k30.meryl
    meryl count compress k=30 threads=XX memory=YY    child.*fastq.gz output    child_compress.k30.meryl

replacing XX and YY with the threads and memory you want meryl to use. Once you have the databases, run:

    $MERQURY/trio/hapmers.sh \
      maternal_compress.k30.meryl \
      paternal_compress.k30.meryl \
         child_compress.k30.meryl

Make sure to count k-mers in compressed space. Child data is optional, in this case, exclude `child_compress.k30.meryl` from the input to `hapmers.sh` and use its output `maternal_compress.k30.only.meryl` and `paternal_compress.k30.only.meryl` in the verkko command below.
</details>

    verkko -d asm \
      --hifi hifi/*.fastq.gz \
      --nano  ont/*.fastq.gz \
      --hap-kmers paternal_compress.k30.hapmer.meryl \
                  maternal_compress.k30.hapmer.meryl \
                  trio


To run in Hi-C mode, reads should be provided using the --hic1 and --hic2 options. For example:

    verkko -d asm \
      --hifi hifi/*.fastq.gz \
      --nano ont/*.fastq.gz \
      --hic1 hic/*R1*fastq.gz  \
      --hic2 hic/*R2*fastq.gz

To run in PoreC mode, reads should be provided using the --porec option. For example:

    verkko -d asm \
      --hifi hifi/*.fastq.gz \
      --nano ont/*.fastq.gz \
      --porec porec/*fastq.gz

Hi-C/PoreC integration was tested mostly on human and primate genomes. Please see  --rdna-tangle, --uneven-depth and --haplo-divergence options if you want to assemble something distant from human and/or have uneven coverage. If you encounter issues or have questions about appropriate parameters, please open an [issue](https://github.com/marbl/verkko/issues).

### Scaffolding:
Verkko includes a separate scaffolding module which is used when Hi-C or Pore-C data are provided (this is separate from Rukki's ability to connect some contigs into scaffolds with just trio)
Verkko tries to makes a rough estimate of the gap size based on the assembly graph. When the gap size estimate is smaller than 100K we report the estimated value. For larger gaps or gaps where a size cannot be estimated, we always report 100K N's.

The scaffolding module uses telomere positions in the assembly (detected with `seqtk telo`), so if your species has a different telomeric repeat motif than vertebrates (CCCTAA), provide it with the `--telomere-motif` option.

If available, you can provide the genome of another individual of the same or closely related species with the `--ref` option. It is _not_ reference-based assembly; the reference will be only used as guidance in scaffolding.

Since the scaffolding module relies on the diploid structure of an assembly, it is not compatible with the `--haploid` option; we recommend the [YaHS](https://github.com/c-zhou/yahs) standalone scaffolder for such cases.

Polyploid scaffolding and phasing is not supported yet.

### Consensus for user-provided paths:

If you already have an assembly but want to customize or change how the nodes in the graph are used, you can do so with the `--paths` option. Verkko will generate contigs from a user-provided file with paths through the assembly graph. Paths should be provided in a [GAF](https://github.com/lh3/gfatools/blob/master/doc/rGFA.md#the-graph-alignment-format-gaf) path format. For example: `name >utig4-1<utig4-2 HAPLOTYPE1` with one path per line, where `utig4-1` is a node in your existing assembly graph and is taken is the fwd orientation and `utig4-2` is in reverse complement. This option also requires the original assembly directory (specified as `--assembly <path_to_assembly>`) and input reads. Output is specified as `-d <output_dir>` (<output_dir> should not be equal to <path_to_assembly>).

### Running on a grid:
By default, verkko will run the snakemake workflow and all compute on the local machine. Support for SGE, Slurm LSF, and PBS (untested) can be enabled with options `--grid`. This will run the snakemake workflow on the local machine but submit all compute to the grid. To launch the both the snakemake workflow and compute on the grid, wrap the verkko command in a shell script and submit using your scheduler. If you're using conda, you may need to make the conda-installed python your default. You can do this with the `--python` option when calling verkko

<details>
<summary>Customizing grid requests (QOS, partition, etc)</summary>
Verkko will submit jobs to the default queue on your grid environment. It is possible to customize how jobs are submitted to specify partitions or other options like accounting or QOS. For example:

```
--snakeopts '--cluster "./slurm-sge-submit.sh {threads} {resources.mem_gb} {resources.time_h} {rulename} {resources.job_id} --partition=quick --account=verkko_asm --qos=verkko_qos"'
```

on SLURM will request the 'quick' queue and pass account and qos options.
</details>

Verkko uses default cpu/memory/time options for different parts of the pipeline. Usually a user does not need to change them, however advanced tuning is possible with `--<stage_code>-run` options.
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

The command line format of all these options is the same: number-of-cpus memory-in-gb time-in-hours, i.e.   `--cns-run 8 32 2`.
Default values can be found in verkko bash script, i.e.  `grep par_ bin/verkko`. Values used for each verkko run are listed in `verkko.yml` in the run directory.
</details>

You can pass through snakemake options to restrict CPU/memory/cluster resources by adding the `--snakeopts` option to verkko. You can also pass arbitrary Snakemake options this way. For example, `--snakeopts "--dry-run"` will print a plan of what jobs will run while `--snakeopts "--cores 1000"` would restrict grid runs to at most 1000 cores across all submited jobs. Running with `--snakeopts "--touch"` will reset all timestamps to be consistent and avoid re-generating any files which already exist. We recommend using `--snakeopts "--dry-run"`, especially when restarting or modifying intermediates in an assembly to check that only the expected steps will run.

### Filtering common contaminants:

Verkko has the ability to filter common contaminants from an assembly using the `--screen` option. You can specify an arbitrary number of targets using multiple `--screen <contaminant_N_name> <contaminant_N_sequence.fasta>`  commands. For each contaminant, verkko will remove all sequences matching the target from the main assembly output. It will also identify a ‘canonical’ representative by coverage and circularize it to remove self-similarity at the start/end. For typical contaminants of human assemblies we have special option `--screen-human-contaminants` (`--screen human` in version 2.2 and earlier) which requires no parameters and is a shortcut for `--screen rDNA rdna.fasta --screen mitochondria mito.fasta --screen EBV ebv.fasta`.
## Outputs:
The final assembly result is under `asm/assembly.fasta`. The final graph (in homopolymer-compressed space) is under `asm/assembly.homopolymer-compressed.gfa` along with coverage files in `asm/assembly*csv`. There is also an `asm/assembly.scfmap` file which translates the final sequence name in `assembly.fasta` to graph nodes. You can find intermediate graphs and coverage files under `asm/*/unitig-*gfa` and `asm/*/unitig-*csv`.

If you provided phasing information, you will also have `asm/assembly.haplotype[12].fasta`, `asm/assembly.colors.csv`, and `asm/assembly.paths.tsv`. The latter two files provide information on the colors obtained from phasing information for each node in the graph and the paths selected to phase the assembly.

If you provided screening information, you will also have an `asm/assembly.exampleN.fasta` and `asm/assembly.exampleN.exemplar.fasta` files which specify all sequences removed from the assembly maching contaminant 'exampleN' as well as the circularized (when possible) selected cannonical sequence.

## Test data:
To test your installation we have an E. coli K12 dataset available.

    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o hifi.fastq.gz
    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o ont.fastq.gz
    verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz

## Citations:
 - Rautiainen M, Nurk S, Walenz BP, Logsdon GA, Porubsky D, Rhie A, Eichler EE, Phillippy AM, Koren S. [Telomere-to-telomere assembly of diploid chromosomes with Verkko](https://doi.org/10.1038/s41587-023-01662-6). Nat Biotech. (2023). `doi:10.1038/s41587-023-01662-6`
 - Antipov D, Rautiainen M, Nurk S, Walenz BP, Solar SJ, Phillippy AM, Koren S. Verkko2 integrates proximity ligation data with long-read De Bruijn graphs for efficient telomere-to-telomere genome assembly, phasing, and scaffolding. Genome Research (2025). `10.1101/gr.280383.124` 
