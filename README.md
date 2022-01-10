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

* Installing with a 'package manager' is encouraged:
  * `conda install -c conda-forge -c bioconda -c defaults MBG`
  * `conda install -c conda-forge -c bioconda -c defaults GraphAligner`
  * `conda install -c conda-forge -c bioconda -c defaults canu`
  * (TODO) Include Canu and Seqrequester as submodules.

* Environment variable VERKKO should be set to the installation directory of Verkko.  This is necessary only if Verkko cannot find it's component scripts.

## Learn:

(To Be Documented)

## Run:

Verkko is implemented as a Snakemake workflow, launched by a wrapper script to parse options
and create a config.yml file.

    verkko -d <work-directory> --hifi <hifi-read-files> --nano <ont-read-files>

Support for SGE, Slurm and LSF (untested) can be enabled with options `--sge`, `--slurm` and `--lsf`, respectively.

## Citations:
 
(In Preparation)
