# Verkko

Verkko is a hybrid genome assembly pipeline developed for
telomere-to-telomere assembly of PacBio HiFi and Oxford Nanopore reads.
Verkko is Finnish for net, mesh and graph.

Verkko uses Canu to correct remaining errors in the HiFi reads, builds a
multiplex de Bruijn graph using MBG, aligns the Oxford Nanopore reads to the
graph using GraphAligner, progressively resolves loops and tangles first with
the HiFi reads then with the aligned Oxford Nanopore reads, and finally
creates contig consensus sequences using Canu's consensus module.

## Install:

* Do NOT download the .zip source code.  It is missing files and will not compile.  This is a [known flaw](https://github.com/dear-github/dear-github/issues/214) with git itself.

* Verkko requires both https://github.com/maickrau/MBG and https://github.com/maickrau/GraphAligner to be installed. Installing with a 'package manager' is encouraged:
  * `conda install -c conda-forge -c bioconda -c defaults MBG`
  * `conda install -c conda-forge -c bioconda -c defaults GraphAligner`

To install Verkko run:

    git clone https://github.com/maickrau/verkko.git
    cd verkko/src
    git submodule init && git submodule update
    make -j32

This will create the folder `verkko/bin` and `verkko/lib/verkko`. You can move the contents of these folders to a central installation location or you can add `verkko/bin` to your path. You may need to set the environment variable VERKKO to the installation directory of Verkko.  This is necessary only if Verkko gives errors that it cannot find component scripts. If GraphAligner or MBG are not available in your path you may also symlink them under `verkko/lib/verkko/bin/`.

## Run:

Verkko is implemented as a Snakemake workflow, launched by a wrapper script to parse options
and create a config.yml file.

    verkko -d <work-directory> --hifi <hifi-read-files> --nano <ont-read-files>

Support for SGE, Slurm and LSF (untested) can be enabled with options `--sge`, `--slurm` and `--lsf`, respectively. To launch on the grid, wrap the verkko command in a shell script and submit to your scheduler.

To test your installation we have an E. coli K12 dataset available. 

    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_hifi_subset24x.fastq.gz -o hifi.fastq.gz
    curl -L https://obj.umiacs.umd.edu/sergek/shared/ecoli_ont_subset50x.fastq.gz -o ont.fastq.gz
    verkko -d asm --hifi ./hifi.fastq.gz --nano ./ont.fastq.gz

The final assembly result is under `asm/7-consensus/concat.fa`. The final graph is under `asm/5-untip/unitig-popped-unitig-normal-connected-tip.gfa` along with coverage files in `asm/5-untip/unitig*csv`. You can find intermediate graphs and coverage files under `asm/*/unitig-*gfa` and `asm/*/unitig-*csv`.

## Citations:
 
(In Preparation)
