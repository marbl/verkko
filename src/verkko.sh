#!/bin/sh

#############################################################################
 #
 #  This file is part of Verkko, a software program that assembles
 #  whole-genome sequencing reads into telomere-to-telomere
 #  haplotype-resolved chromosomes.
 #
 #  Except as indicated otherwise, this is a 'United States Government
 #  Work', and is released in the public domain.
 #
 #  File 'README.licenses' in the root directory of this distribution
 #  contains full conditions and disclaimers.
 #
 ##


#  Use 'cd' and 'pwd' to find the full path to a file/directory.
fullpath() {
  local ipath=$1
  local iname=$(basename $ipath)
  local fpath=

  if [ ! -e $ipath ] ; then
    echo 1>&2 "File '$ipath' doesn't exist."
    exit 1
  fi

  if [ -d $ipath ] ; then      #  If a directory, go in there and
    fpath=$(cd $ipath ; pwd)   #  get the path according to the shell.
  else
    ipath=$(dirname $ipath)    #  If not, strip off the filename,
    fpath=$(cd $ipath ; pwd)   #  then go in, then add the filename
    fpath="$fpath/$iname"      #  back.
  fi

  echo "$fpath"
}



version=""
help=""

hifi=""
nano=""
outd=""
hic1=""
hic2=""

withont="False"
withhic="False"

keepinter="True"

cnspaths=""
cnsdir=""

errors=""
verkko=""

mbg=""
graphaligner=""
mashmap=""
winnowmap=""
python=
perl=

grid="local"
local_cpus=all
local_mem=64

snakeopts=""


#  Set shell variable 'verkko' and environment variable VERKKO to the
#  location of the Verkko components (e.g., /software/verkko/lib/verkko or
#  /usr/local/lib/verkko/).
#
#  If environment variable VERKKO is set, assume that is the path.
#   - this is used when the verkko driver (either src/verkko.sh or
#     bin/verkko) is submitted directly to the grid: the grid will
#     (typically) copy the supplied shell script to a spool directory and run
#     it there - so any auto-detected path based on that path name will be
#     incorrect.
#
#  Otherwise, use the path to this script as a base and look for the
#  lib/verkko components, then set VERKKO in the environment.
#   - if ${verkko}/verkko.sh exists, then we are running out of the source
#     directory (since the 'compiled' version is named 'verkko') and
#     we symlink the compilied binaries here.  This mode of operation
#     will use src/verkko.sh, build/bin/* and src/scripts/* for execution.
#   - if ${verkko}/verkko.sh does not exist, we're running from a fully
#     installed Verkko, and components will be in ${verkko}/../lib/verkko/.
#
#
if [ -e "${VERKKO}/verkko.sh" ] ; then
  verkko=$VERKKO
else
  verkko=$( fullpath $0 )
  verkko=$( dirname $verkko )

  if   [ -e "${verkko}/verkko.sh" ] ; then
    if [ ! -e "${verkko}/bin" ] ; then  cd ${verkko} && ln -s ../build/bin .;  fi
  else
    verkko=$( dirname $verkko )
    verkko="$verkko/lib/verkko"
  fi

  export VERKKO=$verkko
fi

#
#  Algorithm parameters.
#

#  buildStore, countKmers and computeOverlaps
correction_enabled=True

mer_size=28
mer_threshold=20
mer_filter=0.995

cor_min_read=4000
cor_min_overlap=2000
cor_hash_bits=25
cor_filter_kmers=0

#  buildGraph, parameters for MBG
mbg_baseK=1001
mbg_maxK=15000
mbg_window=100
mbg_max_resolution=4000
mbg_hifi_coverage=25
mbg_unitig_abundance=2

#  split_ont, partitioning ont reads for alignment
spl_bases=3000000000
spl_reads=150000
spl_min_length=0

#  align_ont, alignment of ont reads to the initial graph
ali_mxm_length=30
ali_mem_count=10000
ali_bandwidth=15
ali_multi_score_f=0.99
ali_clipping=0.85
ali_min_score=5000
ali_end_clipping=100
ali_incompat_cutoff=0.15
ali_max_trace=5
ali_seed_window=5000

#  post-processing
short_contig_length=100000
screen=""

#  process_ont_paths
pop_min_allowed_cov=5
pop_resolve_steps="20 10 5"
is_haploid="False"

#  Rukki.
ruk_enable="False"
ruk_hap1=""
ruk_hap2=""
ruk_type=""
ruk_fract="0.9"

#  HiC heuristics
haplo_divergence=0.05
#possibly this should be used not only for hic
uneven_depth="False"
no_rdna_tangle="False"



#  split_hic, partitioning illumina hic reads for alignment
#Do not want to use shc_bases at all because of synchronization between left and right that will be broken in case of some precorrecion or trimming
#relying on shc_reads only

shc_bases=30000000000000
shc_reads=20000000
shc_min_length=0



#
#  Run parameters.
#

#  buildStore, countKmers and computeOverlaps
sto_n_cpus=1
sto_mem_gb=4
sto_time_h=4

mer_n_cpus=4
mer_mem_gb=32
mer_time_h=8

meg_n_cpus=4     #  This is temporary, until meryl2 arrives.
meg_mem_gb=32    #  It is used for merging kmer databases.
meg_time_h=4

ovb_n_cpus=8
ovb_mem_gb=32
ovb_time_h=24

ovs_n_cpus=1
ovs_mem_gb=128
ovs_time_h=12

red_n_cpus=4
red_mem_gb=32
red_time_h=4

#  build-graph
mbg_n_cpus=4
mbg_mem_gb=0
mbg_time_h=72

#  process_graph
utg_n_cpus=1
utg_mem_gb=64
utg_time_h=24

#  split_ont
spl_n_cpus=1
spl_mem_gb=8
spl_time_h=96

#  align_ont
ali_n_cpus=24
ali_mem_gb=64
ali_time_h=24

#  process_ont_paths
pop_n_cpus=16
pop_mem_gb=64
pop_time_h=48

#  untip
utp_n_cpus=1
utp_mem_gb=64
utp_time_h=24

#  rukki
ruk_n_cpus=8
ruk_mem_gb=16
ruk_time_h=4

#  create_layout
lay_n_cpus=1
lay_mem_gb=32
lay_time_h=24

#  get_ont_subset
sub_n_cpus=8
sub_mem_gb=16
sub_time_h=24

#  partition_consensus
par_n_cpus=8
par_mem_gb=0
par_time_h=24

#  cns
cns_n_cpus=8
cns_mem_gb=0
cns_time_h=24

# align_hic stuff
ahc_n_cpus=24
ahc_mem_gb=64
ahc_time_h=48

# fast things in hic pipeline
fhc_n_cpus=8
fhc_mem_gb=16
fhc_time_h=24



#
#  If an empty command, give usage.
#

while [ $# -gt 0 ] ; do
    opt=$1   #  Copy 1st word to the 'option' we're testing against.
    arg=$2   #  Copy 2nd word to the (potential) 'argument' for that option.
    arh=$3   #  Copy 3rd word to the (potential) second arg for that option.

    shift    #  Shift out the 'option' word; we're processing it now.

    #
    #  Handle --help and --version special.
    #

    if   [ "$opt" = "--version" ] ;      then echo "$version"; exit 0;
    elif [ "$opt" = "--help" ] ;         then help="help";

    #
    #  Run options
    #

    elif [ "$opt" = "-d" ] ;                   then outd=$arg;          shift
    elif [ "$opt" = "--keep-intermediate" ] ;  then keepinter="True";
    elif [ "$opt" = "--no-cleanup" ] ;         then keepinter="True";
    elif [ "$opt" = "--cleanup" ] ;            then keepinter="False";
    elif [ "$opt" = "--python" ] ;             then python=$arg;        shift
    elif [ "$opt" = "--perl" ] ;               then perl=$arg;          shift
    elif [ "$opt" = "--mbg" ] ;                then mbg=$arg;           shift
    elif [ "$opt" = "--graphaligner" ] ;       then graphaligner=$arg;  shift
    elif [ "$opt" = "--mashmap" ] ;            then mashmap=$arg;       shift
    elif [ "$opt" = "--winnowmap" ] ;          then winnowmap$arg;      shift
    elif [ "$opt" = "--local" ] ;              then grid="local";
    elif [ "$opt" = "--sge" ] ;                then grid="slurm-sge";
    elif [ "$opt" = "--slurm" ] ;              then grid="slurm-sge";
    elif [ "$opt" = "--lsf" ] ;                then grid="lsf";
    elif [ "$opt" = "--local-memory" ] ;       then local_mem=$arg;     shift
    elif [ "$opt" = "--local-cpus" ] ;         then local_cpus=$arg;    shift
    elif [ "$opt" = "--snakeopts" ] ;          then snakeopts=$arg;     shift

    #
    #  Run options for running Rukki before consensus.
    #

    elif [ "$opt" = "--hap-kmers" ] ; then
      ruk_enable="True"
      ruk_hap1=$(fullpath $arg);   shift
      ruk_hap2=$(fullpath $1);     shift
      ruk_type=$1;                 shift

    #
    #  Run options for running only consensus on a set of user-supplied paths.
    #

    elif [ "$opt" = "--paths" ] ;              then cnspaths=$arg;      shift
    elif [ "$opt" = "--assembly" ] ;           then cnsassembly=$arg;   shift

    #
    #  Inputs.
    #
    #  While $arg is a file, append it to the list of input reads.  The test
    #  is checking if the path is absolute or relative; if absolute, the
    #  tested path "//path/to/file" will still point to a file, but if
    #  relative, the tested path "/../subdir/file" will not exist.
    #
    #  It's not a fool-proof method; it is possible that "/data/file" and
    #  "{pwd}/data/file" both exist.
    #
    elif [ "$opt" = "--hifi" ] ; then
        while [ -e "$arg" ] ; do
            if   [ -e "/$arg" ] ; then
                hifi="$hifi $arg"
            else
                hifi="$hifi `pwd`/$arg"
            fi
            shift
            arg=$1
        done

    elif [ "$opt" = "--nano" ] ; then
        withont="True"
        while [ -e "$arg" ] ; do
            if [ -e "/$arg" ] ; then
                nano="$nano $arg"
            else
                nano="$nano `pwd`/$arg"
            fi
            shift
            arg=$1
        done
    elif [ "$opt" = "--hic1" ] ; then
        withhic="True"
        while [ -e "$arg" ] ; do
            if [ -e "/$arg" ] ; then
                hic1="$hic1 $arg"
            else
                hic1="$hic1 `pwd`/$arg"
            fi
            shift
            arg=$1
        done
    elif [ "$opt" = "--hic2" ] ; then
        withhic="True"
        while [ -e "$arg" ] ; do
            if [ -e "/$arg" ] ; then
                hic2="$hic2 $arg"
            else
                hic2="$hic2 `pwd`/$arg"
            fi
            shift
            arg=$1
        done
    elif [ "$opt" = "--no-nano" ] ; then
        withont="False"

    #
    #  Canu correction options
    #

    elif [ "$opt" = "--no-correction" ] ;              then correction_enabled=False

    elif [ "$opt" = "--filter-kmer" ];                 then cor_filter_kmers=$arg; shift
    elif [ "$opt" = "--correct-k-mer-size" ] ;         then mer_size=$arg;        shift
    elif [ "$opt" = "--correct-mer-threshold" ] ;      then mer_threshold=$arg;   shift
    elif [ "$opt" = "--correct-mer-filter" ] ;         then mer_filter=$arg;      shift
    elif [ "$opt" = "--correct-min-read-length" ] ;    then cor_min_read=$arg;    shift
    elif [ "$opt" = "--correct-min-overlap-length" ] ; then cor_min_overlap=$arg; shift
    elif [ "$opt" = "--correct-hash-bits" ] ;          then cor_hash_bits=$arg;   shift

    #
    #  MBG options
    #

    elif [ "$opt" = "--base-k" ] ;        then mbg_baseK=$arg;   shift
    elif [ "$opt" = "--max-k" ] ;         then mbg_maxK=$arg;    shift
    elif [ "$opt" = "--window" ] ;        then mbg_window=$arg;  shift
    elif [ "$opt" = "--threads" ] ;       then mbg_threads=$arg; shift
    elif [ "$opt" = "--max-r" ]   ;       then mbg_max_resolution=$arg; shift
    elif [ "$opt" = "--hifi-coverage" ] ; then mbg_hifi_coverage=$arg; shift
    elif [ "$opt" = "--unitig-abundance" ]; then mbg_unitig_abundance=$arg; shift

    #
    #  splitONT options
    #

    elif [ "$opt" = "--split-bases" ] ;    then spl_bases=$arg;      shift
    elif [ "$opt" = "--split-reads" ] ;    then spl_reads=$arg;      shift
    elif [ "$opt" = "--min-ont-length" ] ; then spl_min_length=$arg; shift

    #
    #  alignONT options
    #

    elif [ "$opt" = "--seed-windows" ];         then ali_seed_window=$arg; shift
    elif [ "$opt" = "--seed-min-length" ] ;     then ali_mxm_length=$arg;      shift
    elif [ "$opt" = "--seed-max-length" ] ;     then ali_mem_count=$arg;       shift
    elif [ "$opt" = "--align-bandwidth" ] ;     then ali_bandwidth=$arg;       shift
    elif [ "$opt" = "--score-fraction" ] ;      then ali_multi_score_f=$arg;   shift
    elif [ "$opt" = "--min-identity" ] ;        then ali_clipping=$arg;        shift
    elif [ "$opt" = "--min-score" ] ;           then ali_min_score=$arg;       shift
    elif [ "$opt" = "--end-clipping" ] ;        then ali_end_clipping=$arg;    shift
    elif [ "$opt" = "--incompatible-cutoff" ] ; then ali_incompat_cutoff=$arg; shift
    elif [ "$opt" = "--max-traces" ] ;          then ali_max_trace=$arg;       shift
    elif [ "$opt" = "--haploid" ];              then is_haploid="True";

    #
    #  HiC options
    #

    elif [ "$opt" = "--no-rdna-tangle" ];          then no_rdna_tangle="True";
    elif [ "$opt" = "--uneven-depth" ];         then uneven_depth="True";
    elif [ "$opt" = "--haplo-divergence" ];     then haplo_divergence=$arg;     shift

    #
    #  Post-processing options
    #

    elif [ "$opt" = "--discard-short" ] ;       then short_contig_length=$arg;   shift
    elif [ "$opt" = "--screen" ] ; then
      if [ "x$arg" = "xhuman" ] ; then
          screen="$screen ebv  human-ebv-AJ507799.2.fasta.gz"
          screen="$screen mito human-mito-NC_012920.1.fasta.gz"
          screen="$screen rdna human-rdna-KY962518.1.fasta.gz"
          shift
      else
          screen="$screen $arg $arh"
          shift   #  Both arg and arh are processed.
          shift
      fi

    #
    #  run-time options
    #

    elif [ "$opt" = "--sto-run" ] ;  then sto_n_cpus=$1; sto_mem_gb=$2; sto_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--mer-run" ] ;  then mer_n_cpus=$1; mer_mem_gb=$2; mer_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--meg-run" ] ;  then meg_n_cpus=$1; meg_mem_gb=$2; meg_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--ovb-run" ] ;  then ovb_n_cpus=$1; ovb_mem_gb=$2; ovb_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--ovs-run" ] ;  then ovs_n_cpus=$1; ovs_mem_gb=$2; ovs_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--red-run" ] ;  then red_n_cpus=$1; red_mem_gb=$2; red_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--mbg-run" ] ;  then mbg_n_cpus=$1; mbg_mem_gb=$2; mbg_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--utg-run" ] ;  then utg_n_cpus=$1; utg_mem_gb=$2; utg_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--spl-run" ] ;  then spl_n_cpus=$1; spl_mem_gb=$2; spl_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--ali-run" ] ;  then ali_n_cpus=$1; ali_mem_gb=$2; ali_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--pop-run" ] ;  then pop_n_cpus=$1; pop_mem_gb=$2; pop_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--utp-run" ] ;  then utp_n_cpus=$1; utp_mem_gb=$2; utp_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--ruk-run" ] ;  then ruk_n_cpus=$1; ruk_mem_gb=$2; ruk_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--lay-run" ] ;  then lay_n_cpus=$1; lay_mem_gb=$2; lay_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--sub-run" ] ;  then sub_n_cpus=$1; sub_mem_gb=$2; sub_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--par-run" ] ;  then par_n_cpus=$1; par_mem_gb=$2; par_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--cns-run" ] ;  then cns_n_cpus=$1; cns_mem_gb=$2; cns_time_h=$3; shift; shift; shift;

    #
    #  unknown options
    #

    else
        errors="${errors}Unknown option '$opt'.\n"
    fi
done

#
#  Set stuff not set by a user-supplied option.
#

if [ "x$mbg" = "x" ] ; then    #  Not set by an option,
  mbg=${verkko}/bin/MBG        #  Set it to our bin/ directory.
fi
if [ ! -e $mbg ] ; then        #  Not in the bin directory,
  mbg=$(which MBG)             #  Set it to whatever is in the PATH.
fi
if [ "x$mbg" != "x" ]; then
  mbg=$(fullpath $mbg)
fi

if [ "x$graphaligner" = "x" ] ; then
  graphaligner=${verkko}/bin/GraphAligner
fi
if [ ! -e $graphaligner ] ; then
  graphaligner=$(which GraphAligner)
fi
if [ "x$graphaligner" != "x" ]; then
  graphaligner=$(fullpath $graphaligner)
fi

if [ "x$mashmap" = "x" ] ; then
  mashmap=${verkko}/bin/mashmap
fi
if [ ! -e $mashmap ] ; then
  mashmap=$(which mashmap)
fi
if [ "x$mashmap" != "x" ]; then
  mashmap=$(fullpath $mashmap)
fi

if [ "x$winnowmap" = "x" ] ; then
  winnowmap=${verkko}/bin/winnowmap
fi
if [ ! -e $winnowmap ] ; then
  winnowmap=$(which winnowmap)
fi
if [ "x$winnowmap" != "x" ]; then
  winnowmap=$(fullpath $winnowmap)
fi
#
#  Fix stuff.
#

if [ "x$cnspaths" != "x" ] ; then
  correction_enabled=False
fi

#
#  Check stuff.
#

if [ "x$outd" = "x" ] ; then
    errors="${errors}No output directory (-d) set.\n"
fi
if [ "x$verkko" = "x" ] ; then
    errors="${errors}Environment variable VERKKO not set.\n"
fi

if [ "x$hifi" = "x" ] ; then
    errors="${errors}No PacBio HiFi reads (--hifi) supplied.\n"
fi

if [ "x$nano" = "x" -a "x$withont" = "xTrue" ] ; then
    errors="${errors}No Oxford Nanopore reads (--nano) supplied.\n"
fi

if [ "x$ruk_enable" = "xTrue" ] ; then
    if [ "x$ruk_hap1" = "x" -o "x$ruk_hap2" = "x" ] ; then
        errors="${errors}Invalid haplotype databases specified, make sure the paths are valid.\n"
    fi
    if [ "x$ruk_type" != "xtrio" -a "x$ruk_type" != "xhic" -a "x$ruk_type" != "xstrandseq" ] ; then
        errors="${errors}Invalid rukki phasing '$ruk_type', must be one of trio/hic/strandseq.\n"
    fi
fi

if [ "x$ruk_enable" = "xTrue" -a -e "$verkko/bin/meryl" ] ; then
    if [ "x$ruk_hap1}" != "x" ] ; then
        $verkko/bin/meryl print threads=1 $ruk_hap1 2> /dev/null | head | grep -q AA
        if [ $? -ne 1 ] ; then
            errors="${errors}Meryl database '$ruk_hap1' appears to be built using non-homopolymer compressed kmers.\n"
        fi
    fi
    if [ "x$ruk_hap2}" != "x" ] ; then
        $verkko/bin/meryl print threads=1 $ruk_hap2 2> /dev/null | head | grep -q AA
        if [ $? -ne 1 ] ; then
            errors="${errors}Meryl database '$ruk_hap2' appears to be built using non-homopolymer compressed kmers.\n"
        fi
    fi
fi

#
#  Check that screen files are present as supplied or in the verkko data directory.
#

if [ ! -z "$screen" ] ; then
    sctest=$screen
    screen=""

    while [ ! -z "$sctest" ]
    do
        sident=$( echo $sctest | cut -s -d ' ' -f 1  )
        sfpath=$( echo $sctest | cut -s -d ' ' -f 2  )
        dtpath="$verkko/data/$sfpath"
        sctest=$( echo $sctest | cut -s -d ' ' -f 3- )

        if   [ -e "$sfpath" ] ; then
            screen="$screen $sident $sfpath"
        elif [ -e "$dtpath" ] ; then
            screen="$screen $sident $dtpath"
        else
            echo "ERROR: Can't find screen '$sident' data file.  Looked for:"
            echo "  user-supplied:   '$sfpath' and"
            echo "  verkko-supplied: '$dtpath'"
        fi
    done
fi

#
#  Check that binaries are present.
#

for exe in bin/findErrors \
           bin/fixErrors \
           bin/layoutToPackage \
           bin/meryl \
           bin/meryl-lookup \
           bin/ovStoreBuild \
           bin/ovStoreConfig \
           bin/overlapInCore \
           bin/overlapInCorePartition \
           bin/sqStoreCreate \
           bin/sqStoreDumpMetaData \
           bin/utgcns \
           scripts/get_layout_from_mbg.py ; do
  if [ ! -e "$verkko/$exe" ] ; then
      errors="${errors}Can't find '$exe' in directory VERKKO = '$verkko/'.\n"
  fi
done

if   [ "x$mbg" = "x" ] ; then
    errors="${errors}Can't find MBG executable in \$PATH or \$VERKKO/bin/MBG.\n"
elif [ ! -e "$mbg" ] ; then
    errors="${errors}Can't find MBG executable at '$mbg'.\n"
fi

if   [ "x$graphaligner" = "x" ] ; then
    errors="${errors}Can't find GraphAligner executable in \$PATH or \$VERKKO/bin/GraphAligner.\n"
elif [ ! -e "$graphaligner" ] ; then
    errors="${errors}Can't find GraphAligner executable at '$graphaligner'.\n"
fi

if   [ "x$mashmap" = "x" ] ; then
    errors="${errors}Can't find mashmap executable in \$PATH or \$VERKKO/bin/mashmap.\n"
elif [ ! -e "$mashmap" ] ; then
    errors="${errors}Can't find mashmap executable at '$mashmap'.\n"
fi

if   [ "x$winnowmap" = "x" ] ; then
    errors="${errors}Can't find winnowmap executable in \$PATH or \$VERKKO/bin/winnowmap.\n"
elif [ ! -e "$winnowmap" ] ; then
    errors="${errors}Can't find winnowmap executable at '$winnowmap'.\n"
fi

#
#  Complain!
#

if [ "x$help" = "xhelp" -o "x$errors" != "x" ] ; then
    echo "usage: $0 -d <output-directory> --hifi <hifi-reads ...> --nano <nanopore-reads ...>"
    echo "  MANDATORY PARAMETERS:"
    echo "    -d <output-directory>    Directory to use for verkko intermediate and final results."
    echo "                             Will be created if needed."
    echo "    --hifi <files ...>       List of files containing PacBio HiFi reads."
    echo "    --nano <files ...>       List of files containing Oxford Nanopore reads."
    echo ""
    echo "                             Input reads can be any combination of FASTA/FASTQ,"
    echo "                             uncompressed or gzip/bzip2/xz compressed.  Any"
    echo "                             number of files can be supplied; *.gz works."
    echo ""
    echo "  ALGORITHM PARAMETERS:"
    echo "    --no-correction          Do not perform Canu correction on the HiFi reads."
    echo "    --no-nano                Assemble without ONT data."
    echo ""
    echo "    --hap-kmers h1 h2 type   Use rukki to assign paths to haplotypes.  'h1' and 'h2"
    echo "                             must be Meryl databases of homopolymer-compressed parental"
    echo "                             kmers.  'type' must be 'trio', 'hic' or 'strandseq'."
    echo ""
    echo "    --hic1 <files ...>       List of files containing left hic reads."
    echo "    --hic2 <files ...>       List of files containing right hic reads."
    echo "                             Order of left and right files should be same, no interlaced files allowed."
    echo "                             Input reads can be any combination of FASTA/FASTQ,"
    echo "                             uncompressed or gzip/bzip2/xz compressed.  Any"
    echo "                             number of files can be supplied; *.gz works."
    echo "    --no-rdna-tangle         Switch off option that helps to proceed large rDNA tangles which may connect multiple chromosomes."
    echo "    --uneven-depth           Disable coverage-based heuristics in homozygous nodes detection for phasing."
    echo "    --haplo-divergence       Estimation on maximum divergence between haplotypes. Should be increased for species with divergence significantly higher than in human. Default: 0.05, min 0, max 0.2"
    echo ""
    echo "    --screen <option>        Identify common contaminants and remove from the assembly, saving 1 (circularized) exemplar."
    echo "                             For human, '--screen human' will attempt to remove rDNA, mitochondria, and EBV."
    echo "                             Arbitrary contaminants are supported by supplying a name and fasta:"
    echo "                             'screen contaminant /full/path/to/contaminant.fasta'"
    echo "                             Multiple screen commands are allowed and are additive."
    echo ""
    echo "    --paths <gaf paths>      No assembly, generate consensus given paths and an existing assembly."
    echo "                             The gaf file must be formatted as follows: 'name >utig4-1<utig4-2 HAPLOTYPE1'. One per line."
    echo "                             where utig4-1 is in fwd orientation and utig4-2 is in reverse complement. Requires '--assembly'."
    echo "    --assembly <output-dir>  Existing verkko assembly where the nodes given to --paths are defined."
    echo "                             The nodes should come from assembly.homopolymer-compressed.gfa"
    echo "                             Provide -d for the output of new consensus, as well as the hifi and nano reads from previous run."
    echo ""
    echo "  COMPUTATIONAL PARAMETERS:"
    echo "    --python <interpreter>   Path or name of a python interpreter.  Default: 'python'."
    echo "    --perl <interpreter>     Path of name of a perl interpreter.  Default: 'perl'."
    echo "    --mbg <path>             Path to MBG.             Default for all three"
    echo "    --graphaligner <path>    Path to GraphAligner.    one packaged with verkko."
    echo "    --mashmap <path>         Path to mashmap."
    echo "    --winnowmap <path>       Path to winnowmap."
    echo ""
    echo "    --cleanup                Remove intermediate results."
    echo "    --no-cleanup             Retain intermediate results (default)."
    echo ""
    echo "    --local                  Run on the local machine (default)."
    echo "    --local-memory           Specify the upper limit on memory to use, in GB, default 64"
    echo "    --local-cpus             Specify the number of CPUs to use, default 'all'"
    echo ""
    echo "    --sge                    Enable Sun Grid Engine support."
    echo "    --slurm                  Enable Slurm support."
    echo "    --lsf                    Enable IBM Spectrum LSF support."
    echo ""
    echo "    --snakeopts <string>     Append snakemake options in \"string\" to the"
    echo "                             snakemake command.  Options MUST be quoted."
    echo "                             Can be used to change default resources, or"
    echo "                             pass additional parameters to the grid engine"
    echo "                             by setting the resource 'clusterparams'."
    echo ""
    echo "    --sto-run                Set resource limits for various stages."
    echo "    --mer-run                Format: number-of-cpus memory-in-gb time-in-hours"
    echo "    --ovb-run                  --cns-run 8 32 2"
    echo "    --ovs-run"
    echo "    --red-run"
    echo "    --mbg-run"
    echo "    --utg-run"
    echo "    --spl-run"
    echo "    --ali-run"
    echo "    --pop-run"
    echo "    --utp-run"
    echo "    --lay-run"
    echo "    --sub-run"
    echo "    --par-run"
    echo "    --cns-run"
    echo ""
    echo "ADVANCED MODULE PARAMETERS (expert users):"
    echo "HiFi read correction:"
    echo "    --correct-k-mer-size"
    echo "    --correct-mer-threshold"
    echo "    --correct-mer-filter"
    echo "    --correct-min-read-length"
    echo "    --correct-min-overlap-length"
    echo "    --correct-hash-bits"
    echo "    "
    echo "MBG:"
    echo "    --base-k"
    echo "    --max-k"
    echo "    --window"
    echo "    --threads"
    echo "    --unitig-abundance"
    echo "    --hifi-coverage"
    echo "    "
    echo "ONT splitting:"
    echo "    --split-bases"
    echo "    --split-reads"
    echo "    --min-ont-length"
    echo "    "
    echo "    "
    echo "GraphAligner:"
    echo "    --seed-min-length"
    echo "    --seed-max-length"
    echo "    --align-bandwidth"
    echo "    --score-fraction"
    echo "    --min-identity"
    echo "    --min-score"
    echo "    --end-clipping"
    echo "    --incompatible-cutoff"
    echo "    --max-trace"
    echo ""
    echo "  Verkko module path: ${verkko}/"
    echo ""

    if [ "x$help" = "xhelp" ]
    then       #  Requested --help, not an error.
        exit 0
    else       #  Errors!
        printf "${errors}"
        exit 1
    fi
fi

#
#  All good!
#


#
#  Make a YAML configuration for Snakemake.
#

mkdir -p ${outd}

echo  > ${outd}/verkko.yml "#  Generated automatically by verkko.sh."
echo >> ${outd}/verkko.yml "#  Changes will be overwritten."
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "VERKKO:              '${verkko}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "MBG:                 '${mbg}'"
echo >> ${outd}/verkko.yml "GRAPHALIGNER:        '${graphaligner}'"
echo >> ${outd}/verkko.yml "MASHMAP:             '${mashmap}'"
echo >> ${outd}/verkko.yml "WINNOWMAP:           '${winnowmap}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "PYTHON:              '${python}'"
echo >> ${outd}/verkko.yml "PERL:                '${perl}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "HIFI_READS:"
for h in ${hifi} ; do
  echo >> ${outd}/verkko.yml " - '$h'"
done
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "withONT:             '${withont}'"
echo >> ${outd}/verkko.yml "ONT_READS:"
for o in ${nano} ; do
  echo >> ${outd}/verkko.yml " - '$o'"
done
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "withHIC:             '${withhic}'"
echo >> ${outd}/verkko.yml "HIC_READS1:"
for o in ${hic1} ; do
  echo >> ${outd}/verkko.yml " - '$o'"
done
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "HIC_READS2:"
for o in ${hic2} ; do
  echo >> ${outd}/verkko.yml " - '$o'"
done
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  Algorithm parameters."
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  buildStore, countKmers and computeOverlaps"
echo >> ${outd}/verkko.yml "correction_enabled:  '${correction_enabled}'"
echo >> ${outd}/verkko.yml "mer_size:            '${mer_size}'"
echo >> ${outd}/verkko.yml "mer_threshold:       '${mer_threshold}'"
echo >> ${outd}/verkko.yml "mer_filter:          '${mer_filter}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "cor_min_read:        '${cor_min_read}'"
echo >> ${outd}/verkko.yml "cor_min_overlap:     '${cor_min_overlap}'"
echo >> ${outd}/verkko.yml "cor_hash_bits:       '${cor_hash_bits}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "cor_filter_kmers:    '${cor_filter_kmers}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  build-graph, MBG"
echo >> ${outd}/verkko.yml "mbg_baseK:           '${mbg_baseK}'"
echo >> ${outd}/verkko.yml "mbg_maxK:            '${mbg_maxK}'"
echo >> ${outd}/verkko.yml "mbg_window:          '${mbg_window}'"
echo >> ${outd}/verkko.yml "mbg_max_resolution:  '${mbg_max_resolution}'"
echo >> ${outd}/verkko.yml "mbg_hifi_coverage:   '${mbg_hifi_coverage}'"
echo >> ${outd}/verkko.yml "mbg_unitig_abundance: '${mbg_unitig_abundance}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  split_ont"
echo >> ${outd}/verkko.yml "spl_bases:           '${spl_bases}'"
echo >> ${outd}/verkko.yml "spl_reads:           '${spl_reads}'"
echo >> ${outd}/verkko.yml "spl_min_length:      '${spl_min_length}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  align_ont, GraphAligner"
echo >> ${outd}/verkko.yml "ali_mxm_length:      '${ali_mxm_length}'"
echo >> ${outd}/verkko.yml "ali_mem_count:       '${ali_mem_count}'"
echo >> ${outd}/verkko.yml "ali_bandwidth:       '${ali_bandwidth}'"
echo >> ${outd}/verkko.yml "ali_multi_score_f:   '${ali_multi_score_f}'"
echo >> ${outd}/verkko.yml "ali_clipping:        '${ali_clipping}'"
echo >> ${outd}/verkko.yml "ali_min_score:       '${ali_min_score}'"
echo >> ${outd}/verkko.yml "ali_end_clipping:    '${ali_end_clipping}'"
echo >> ${outd}/verkko.yml "ali_incompat_cutoff: '${ali_incompat_cutoff}'"
echo >> ${outd}/verkko.yml "ali_max_trace:       '${ali_max_trace}'"
echo >> ${outd}/verkko.yml "ali_seed_window:     '${ali_seed_window}'"
echo >> ${outd}/verkko.yml "cor_filter_kmers:    '${cor_filter_kmers}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  post-processing"
echo >> ${outd}/verkko.yml "short_contig_length: '${short_contig_length}'"
echo >> ${outd}/verkko.yml "screen:              '${screen}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  process_ont_paths"
echo >> ${outd}/verkko.yml "pop_min_allowed_cov: '${pop_min_allowed_cov}'"
echo >> ${outd}/verkko.yml "pop_resolve_steps:   '${pop_resolve_steps}'"
echo >> ${outd}/verkko.yml "haploid:             '${is_haploid}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  Rukki"
echo >> ${outd}/verkko.yml "ruk_enable:          '${ruk_enable}'"
echo >> ${outd}/verkko.yml "ruk_hap1:            '${ruk_hap1}'"
echo >> ${outd}/verkko.yml "ruk_hap2:            '${ruk_hap2}'"
echo >> ${outd}/verkko.yml "ruk_type:            '${ruk_type}'"
echo >> ${outd}/verkko.yml "ruk_fract:           '${ruk_fract}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  HiC algo options"
echo >> ${outd}/verkko.yml "no_rdna_tangle:      '${no_rdna_tangle}'"
echo >> ${outd}/verkko.yml "uneven_depth:        '${uneven_depth}'"
echo >> ${outd}/verkko.yml "haplo_divergence:    '${haplo_divergence}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  Aligning hic reads"
echo >> ${outd}/verkko.yml "shc_bases:           '${shc_bases}'"
echo >> ${outd}/verkko.yml "shc_reads:           '${shc_reads}'"
echo >> ${outd}/verkko.yml "shc_min_length:      '${shc_min_length}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  Run parameters."
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "keep_intermediate:   '${keepinter}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  buildStore, countKmers and computeOverlaps"
echo >> ${outd}/verkko.yml "sto_n_cpus:          '${sto_n_cpus}'"
echo >> ${outd}/verkko.yml "sto_mem_gb:          '${sto_mem_gb}'"
echo >> ${outd}/verkko.yml "sto_time_h:          '${sto_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "mer_n_cpus:          '${mer_n_cpus}'"
echo >> ${outd}/verkko.yml "mer_mem_gb:          '${mer_mem_gb}'"
echo >> ${outd}/verkko.yml "mer_time_h:          '${mer_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "meg_n_cpus:          '${meg_n_cpus}'"
echo >> ${outd}/verkko.yml "meg_mem_gb:          '${meg_mem_gb}'"
echo >> ${outd}/verkko.yml "meg_time_h:          '${meg_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "ovb_n_cpus:          '${ovb_n_cpus}'"
echo >> ${outd}/verkko.yml "ovb_mem_gb:          '${ovb_mem_gb}'"
echo >> ${outd}/verkko.yml "ovb_time_h:          '${ovb_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "ovs_n_cpus:          '${ovs_n_cpus}'"
echo >> ${outd}/verkko.yml "ovs_mem_gb:          '${ovs_mem_gb}'"
echo >> ${outd}/verkko.yml "ovs_time_h:          '${ovs_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "red_n_cpus:          '${red_n_cpus}'"
echo >> ${outd}/verkko.yml "red_mem_gb:          '${red_mem_gb}'"
echo >> ${outd}/verkko.yml "red_time_h:          '${red_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  build-graph"
echo >> ${outd}/verkko.yml "mbg_n_cpus:          '${mbg_n_cpus}'"
echo >> ${outd}/verkko.yml "mbg_mem_gb:          '${mbg_mem_gb}'"
echo >> ${outd}/verkko.yml "mbg_time_h:          '${mbg_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  process_graph"
echo >> ${outd}/verkko.yml "utg_n_cpus:          '${utg_n_cpus}'"
echo >> ${outd}/verkko.yml "utg_mem_gb:          '${utg_mem_gb}'"
echo >> ${outd}/verkko.yml "utg_time_h:          '${utg_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  split_ont"
echo >> ${outd}/verkko.yml "spl_n_cpus:          '${spl_n_cpus}'"
echo >> ${outd}/verkko.yml "spl_mem_gb:          '${spl_mem_gb}'"
echo >> ${outd}/verkko.yml "spl_time_h:          '${spl_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  align_ont"
echo >> ${outd}/verkko.yml "ali_n_cpus:          '${ali_n_cpus}'"
echo >> ${outd}/verkko.yml "ali_mem_gb:          '${ali_mem_gb}'"
echo >> ${outd}/verkko.yml "ali_time_h:          '${ali_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  process_ont_paths"
echo >> ${outd}/verkko.yml "pop_n_cpus:          '${pop_n_cpus}'"
echo >> ${outd}/verkko.yml "pop_mem_gb:          '${pop_mem_gb}'"
echo >> ${outd}/verkko.yml "pop_time_h:          '${pop_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  untip"
echo >> ${outd}/verkko.yml "utp_n_cpus:          '${utp_n_cpus}'"
echo >> ${outd}/verkko.yml "utp_mem_gb:          '${utp_mem_gb}'"
echo >> ${outd}/verkko.yml "utp_time_h:          '${utp_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  rukki"
echo >> ${outd}/verkko.yml "ruk_n_cpus:          '${ruk_n_cpus}'"
echo >> ${outd}/verkko.yml "ruk_mem_gb:          '${ruk_mem_gb}'"
echo >> ${outd}/verkko.yml "ruk_time_h:          '${ruk_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  create_layout"
echo >> ${outd}/verkko.yml "lay_n_cpus:          '${lay_n_cpus}'"
echo >> ${outd}/verkko.yml "lay_mem_gb:          '${lay_mem_gb}'"
echo >> ${outd}/verkko.yml "lay_time_h:          '${lay_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  get_ont_subset"
echo >> ${outd}/verkko.yml "sub_n_cpus:          '${sub_n_cpus}'"
echo >> ${outd}/verkko.yml "sub_mem_gb:          '${sub_mem_gb}'"
echo >> ${outd}/verkko.yml "sub_time_h:          '${sub_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  partition_consensus"
echo >> ${outd}/verkko.yml "par_n_cpus:          '${par_n_cpus}'"
echo >> ${outd}/verkko.yml "par_mem_gb:          '${par_mem_gb}'"
echo >> ${outd}/verkko.yml "par_time_h:          '${par_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  cns"
echo >> ${outd}/verkko.yml "cns_n_cpus:          '${cns_n_cpus}'"
echo >> ${outd}/verkko.yml "cns_mem_gb:          '${cns_mem_gb}'"
echo >> ${outd}/verkko.yml "cns_time_h:          '${cns_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  align hic stuff"
echo >> ${outd}/verkko.yml "ahc_n_cpus:          '${ahc_n_cpus}'"
echo >> ${outd}/verkko.yml "ahc_mem_gb:          '${ahc_mem_gb}'"
echo >> ${outd}/verkko.yml "ahc_time_h:          '${ahc_time_h}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  fast hic stuff"
echo >> ${outd}/verkko.yml "fhc_n_cpus:          '${fhc_n_cpus}'"
echo >> ${outd}/verkko.yml "fhc_mem_gb:          '${fhc_mem_gb}'"
echo >> ${outd}/verkko.yml "fhc_time_h:          '${fhc_time_h}'"

echo >> ${outd}/verkko.yml "#  This is the end."

#
#  If cnspaths (and cnsassembly) are defined and exist, set up to run just
#  consensus on the paths in ${cnspaths}, using various bits from an existing
#  Verkko assembly in ${cnsassembly}.
#

target="verkko"

if [ "x$withhic" = "xTrue" ] ; then
    target="runRukkiHIC"
fi

if [ "x$cnspaths" != "x" ] ; then
    target="cnspath"

    cnspaths=$(fullpath $cnspaths)         #  Convert to full
    cnsassembly=$(fullpath $cnsassembly)   #  path.

    if [ ! -e $cnspaths ] ; then
        echo "Can't find --paths ${cnspaths}."
        exit 1
    fi

    if [ ! -d $cnsassembly ] ; then
        echo "Can't find --assembly ${cnsassembly}."
        exit 1
    fi

    #  Copy pieces from the previous assembly to the new run directory.  This
    #  is done - instead of symlinking - to prevent Snakemake from
    #  'accidentally' obliterating precious original files.

    if [ ! -e "${outd}/5-untip" ]; then
       mkdir ${outd}/5-untip
       cp -p ${cnsassembly}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv     ${outd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv
       cp -p ${cnsassembly}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.ont-coverage.csv      ${outd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.ont-coverage.csv
       cp -p ${cnsassembly}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa                   ${outd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa
       cp -p ${cnsassembly}/5-untip/combined-edges-final.gfa                                        ${outd}/5-untip/combined-edges-final.gfa
       cp -p ${cnsassembly}/5-untip/combined-nodemap-final.txt                                      ${outd}/5-untip/combined-nodemap-final.txt
       cp -p ${cnsassembly}/5-untip/nodelens-final.txt                                              ${outd}/5-untip/nodelens-final.txt
    fi

    if [ ! -e "${outd}/6-layoutContigs" ] ; then
        mkdir ${outd}/6-layoutContigs
        cp -p ${cnsassembly}/6-layoutContigs/combined-nodemap.txt     ${outd}/6-layoutContigs/combined-nodemap.txt
        cp -p ${cnsassembly}/6-layoutContigs/combined-edges.gfa       ${outd}/6-layoutContigs/combined-edges.gfa
        cp -p ${cnsassembly}/6-layoutContigs/combined-alignments.gaf  ${outd}/6-layoutContigs/combined-alignments.gaf
        cp -p ${cnspaths}                                             ${outd}/6-layoutContigs/consensus_paths.txt
        cp -p ${cnsassembly}/6-layoutContigs/nodelens.txt             ${outd}/6-layoutContigs/nodelens.txt
    fi

    if [ ! -e "${outd}/7-consensus" ] ; then
        mkdir ${outd}/7-consensus
        cp -p ${cnsassembly}/7-consensus/ont_subset.fasta.gz          ${outd}/7-consensus/ont_subset.fasta.gz
    fi
fi

#
#  Generate a script to launch snakemake.
#

#  Detect if the '--rerun-triggers mtime' is valid, and enable it.  This fixes
#  a serious incompatibility in snakemake 7.8.x.
sv=$(snakemake --version)
svmaj=$(echo $sv | cut -d. -f 1)
svmin=$(echo $sv | cut -d. -f 2)
svpat=$(echo $sv | cut -d. -f 3)

echo  > ${outd}/snakemake.sh "#!/bin/sh"
echo >> ${outd}/snakemake.sh ""
echo >> ${outd}/snakemake.sh "echo \"Launching $version\""
echo >> ${outd}/snakemake.sh "echo \"Using snakemake $svmaj.$svmin.$svpat.\""
echo >> ${outd}/snakemake.sh ""
echo >> ${outd}/snakemake.sh "snakemake ${target} --nocolor \\"
echo >> ${outd}/snakemake.sh "  --directory . \\"
echo >> ${outd}/snakemake.sh "  --snakefile ${verkko}/Snakefile \\"
echo >> ${outd}/snakemake.sh "  --configfile verkko.yml \\"
echo >> ${outd}/snakemake.sh "  --reason \\"
echo >> ${outd}/snakemake.sh "  --keep-going \\"
echo >> ${outd}/snakemake.sh "  --rerun-incomplete \\"
if [ $svmaj -ge 7 -a $svmin -ge 8 ] ; then
    echo >> ${outd}/snakemake.sh "  --rerun-triggers mtime \\"
fi
if [ $grid = "local" ] ; then
    echo >> ${outd}/snakemake.sh "  --latency-wait 2 \\"
    echo >> ${outd}/snakemake.sh "  --cores ${local_cpus} \\"
    echo >> ${outd}/snakemake.sh "  --resources mem_gb=${local_mem} \\"
else
    echo >> ${outd}/snakemake.sh "  --latency-wait 30 \\"
    echo >> ${outd}/snakemake.sh "  --jobs 1000 \\"
    echo >> ${outd}/snakemake.sh "  --profile ${verkko}/profiles \\"
    echo >> ${outd}/snakemake.sh "  --restart-times 1 \\"
    echo >> ${outd}/snakemake.sh "  --max-jobs-per-second 10 \\"
    echo >> ${outd}/snakemake.sh "  --max-status-checks-per-second 0.02 \\"
    echo >> ${outd}/snakemake.sh "  --local-cores 1 \\"
fi
echo >> ${outd}/snakemake.sh "  ${snakeopts}"
echo >> ${outd}/snakemake.sh ""

chmod +x ${outd}/snakemake.sh

#
#  Move into the output directory and run the command.
#

cd ${outd}
./snakemake.sh


#Failed to do it with snakemake
if [ "x$withhic" = "xTrue" ] ; then
    if [ ! -e "8-hicPipeline/rukki.paths.gaf" ]; then
        echo "ERROR!"
        echo "Not running final consensus since no rukki paths provided!"
        exit
    fi
    newoutd=8-hicPipeline/final_contigs/
    mkdir -p $newoutd
    cp verkko.yml $newoutd
    cp snakemake.sh $newoutd
    if [ ! -e "${newoutd}/5-untip" ]; then
       mkdir ${newoutd}/5-untip
       cp -p 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv     ${newoutd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.hifi-coverage.csv
       cp -p 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.ont-coverage.csv      ${newoutd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.ont-coverage.csv
       cp -p 5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa                   ${newoutd}/5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa
       cp -p 5-untip/combined-edges-final.gfa                                        ${newoutd}/5-untip/combined-edges-final.gfa
       cp -p 5-untip/combined-nodemap-final.txt                                      ${newoutd}/5-untip/combined-nodemap-final.txt
       cp -p 5-untip/nodelens-final.txt                                              ${newoutd}/5-untip/nodelens-final.txt
    fi

    if [ ! -e "${newoutd}/6-layoutContigs" ] ; then
        mkdir ${newoutd}/6-layoutContigs
        cp -p 6-layoutContigs/combined-nodemap.txt     ${newoutd}/6-layoutContigs/combined-nodemap.txt
        cp -p 6-layoutContigs/combined-edges.gfa       ${newoutd}/6-layoutContigs/combined-edges.gfa
        cp -p 6-layoutContigs/combined-alignments.gaf  ${newoutd}/6-layoutContigs/combined-alignments.gaf
        cp -p 8-hicPipeline/rukki.paths.gaf            ${newoutd}/6-layoutContigs/consensus_paths.txt
        cp -p 6-layoutContigs/nodelens.txt             ${newoutd}/6-layoutContigs/nodelens.txt
    fi

    if [ ! -e "${newoutd}/7-consensus" ] ; then
        mkdir ${newoutd}/7-consensus
        cp -p 7-consensus/ont_subset.fasta.gz          ${newoutd}/7-consensus/ont_subset.fasta.gz
    fi
    cd $newoutd
    sed -i 's/runRukkiHIC/cnspath/g' snakemake.sh
    ./snakemake.sh
    cp *.fasta ../../
    cp *.layout ../../
fi
