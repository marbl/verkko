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
    exit
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

withont="False"

keepinter="False"

cnspaths=""
cnsdir=""

errors=""
verkko=""

mbg=""
graphaligner=""
python=

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

cor_min_read=4000
cor_min_overlap=2000
cor_hash_bits=25

#  buildGraph, parameters for MBG
mbg_baseK=1001
mbg_maxK=15000
mbg_window=100
mbg_max_resolution=4000

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

#  process_ont_paths
pop_min_allowed_cov=5
pop_resolve_steps="20 10 5"

#  Rukki.
ruk_enable="False"
ruk_hap1=""
ruk_hap2=""
ruk_type=""
ruk_fract="0.9"

#
#  Run parameters.
#

#  buildStore, countKmers and computeOverlaps
sto_n_cpus=1
sto_mem_gb=4
sto_time_h=4

mer_n_cpus=4
mer_mem_gb=16
mer_time_h=4

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
mbg_mem_gb=128
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
pop_n_cpus=1
pop_mem_gb=64
pop_time_h=24

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

#
#  If an empty command, give usage.
#

while [ $# -gt 0 ] ; do
    opt=$1
    arg=$2

    shift

    #
    #  Handle --help and --version special.
    #

    if   [ "$opt" = "--version" ] ;      then echo "$version"; exit;
    elif [ "$opt" = "--help" ] ;         then help="help";

    #
    #  Run options
    #

    elif [ "$opt" = "-d" ] ;                   then outd=$arg;          shift
    elif [ "$opt" = "--keep-intermediate" ] ;  then keepinter="True";
    elif [ "$opt" = "--python" ] ;             then python=$arg;        shift
    elif [ "$opt" = "--mbg" ] ;                then mbg=$arg;           shift
    elif [ "$opt" = "--graphaligner" ] ;       then graphaligner=$arg;  shift
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

    elif [ "$opt" = "--no-nano" ] ; then
        withont="False"

    #
    #  Canu correction options
    #

    elif [ "$opt" = "--no-correction" ] ;              then correction_enabled=False

    elif [ "$opt" = "--correct-k-mer-size" ] ;         then mer_size=$arg;        shift
    elif [ "$opt" = "--correct-mer-threshold" ] ;      then mer_threshold=$arg;   shift
    elif [ "$opt" = "--correct-min-read-length" ] ;    then cor_min_read=$arg;    shift
    elif [ "$opt" = "--correct-min-overlap-length" ] ; then cor_min_overlap=$arg; shift
    elif [ "$opt" = "--correct-hash-bits" ] ;          then cor_hash_bits=$arg;   shift 

    #
    #  MBG options
    #

    elif [ "$opt" = "--base-k" ] ;   then mbg_baseK=$arg;   shift
    elif [ "$opt" = "--max-k" ] ;    then mbg_maxK=$arg;    shift
    elif [ "$opt" = "--window" ] ;   then mbg_window=$arg;  shift
    elif [ "$opt" = "--threads" ] ;  then mbg_threads=$arg; shift
    elif [ "$opt" = "--max-r" ]   ;  then mbg_max_resolution=$arg; shift

    #
    #  splitONT options
    #

    elif [ "$opt" = "--split-bases" ] ;    then spl_bases=$arg;      shift
    elif [ "$opt" = "--split-reads" ] ;    then spl_reads=$arg;      shift
    elif [ "$opt" = "--min-ont-length" ] ; then spl_min_length=$arg; shift

    #
    #  alignONT options
    #

    elif [ "$opt" = "--seed-min-length" ] ;     then ali_mxm_length=$arg;      shift
    elif [ "$opt" = "--seed-max-length" ] ;     then ali_mem_count=$arg;       shift
    elif [ "$opt" = "--align-bandwidth" ] ;     then ali_bandwidth=$arg;       shift
    elif [ "$opt" = "--score-fraction" ] ;      then ali_multi_score_f=$arg;   shift
    elif [ "$opt" = "--min-identity" ] ;        then ali_clipping=$arg;        shift
    elif [ "$opt" = "--min-score" ] ;           then ali_min_score=$arg;       shift
    elif [ "$opt" = "--end-clipping" ] ;        then ali_end_clipping=$arg;    shift
    elif [ "$opt" = "--incompatible-cutoff" ] ; then ali_incompat_cutoff=$arg; shift
    elif [ "$opt" = "--max-traces" ] ;          then ali_max_trace=$arg;       shift

    #
    #  run-time options
    #

    elif [ "$opt" = "--sto-run" ] ;  then sto_n_cpus=$1; sto_mem_gb=$2; sto_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--mer-run" ] ;  then mer_n_cpus=$1; mer_mem_gb=$2; mer_time_h=$3; shift; shift; shift;
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
    echo "    --hap-kmers h1 h2 type  Use rukki to assign paths to haplotypes.  'h1' and 'h2"
    echo "                            must be Meryl databases of homopolymer-compressed parental"
    echo "                            kmers.  'type' must be 'trio', 'hic' or 'strandseq'."
    echo ""
    echo "    --base-k"
    echo "    --max-k"
    echo "    --window"
    echo "    --threads"
    echo "    "
    echo "    --split-bases"
    echo "    --split-reads"
    echo "    --min-ont-length"
    echo "    "
    echo "    --correct-k-mer-size"
    echo "    --correct-mer-threshold"
    echo "    --correct-min-read-length"
    echo "    --correct-min-overlap-length"
    echo "    --correct-hash-bits" 
    echo "    "
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
    echo "  COMPUTATIONAL PARAMETERS:"
    echo "    --python <interpreter>   Path or name of a python interpreter.  Default: 'python'."
    echo "    --mbg <path>             Path to MBG.             Default for both is the"
    echo "    --graphaligner <path>    Path to GraphAligner.    one packaged with verkko."
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
    echo "  Verkko module path: ${verkko}/"
    echo ""
    printf "${errors}"
    exit 0
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
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "PYTHON:              '${python}'"
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
echo >> ${outd}/verkko.yml "#  Algorithm parameters."
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  buildStore, countKmers and computeOverlaps"
echo >> ${outd}/verkko.yml "correction_enabled:  '${correction_enabled}'"
echo >> ${outd}/verkko.yml "mer_size:            '${mer_size}'"
echo >> ${outd}/verkko.yml "mer_threshold:       '${mer_threshold}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "cor_min_read:        '${cor_min_read}'"
echo >> ${outd}/verkko.yml "cor_min_overlap:     '${cor_min_overlap}'"
echo >> ${outd}/verkko.yml "cor_hash_bits:       '${cor_hash_bits}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  build-graph, MBG"
echo >> ${outd}/verkko.yml "mbg_baseK:           '${mbg_baseK}'"
echo >> ${outd}/verkko.yml "mbg_maxK:            '${mbg_maxK}'"
echo >> ${outd}/verkko.yml "mbg_window:          '${mbg_window}'"
echo >> ${outd}/verkko.yml "mbg_max_resolution:  '${mbg_max_resolution}'"
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
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  process_ont_paths"
echo >> ${outd}/verkko.yml "pop_min_allowed_cov: '${pop_min_allowed_cov}'"
echo >> ${outd}/verkko.yml "pop_resolve_steps:   '${pop_resolve_steps}'"
echo >> ${outd}/verkko.yml ""
echo >> ${outd}/verkko.yml "#  Rukki"
echo >> ${outd}/verkko.yml "ruk_enable:          '${ruk_enable}'"
echo >> ${outd}/verkko.yml "ruk_hap1:            '${ruk_hap1}'"
echo >> ${outd}/verkko.yml "ruk_hap2:            '${ruk_hap2}'"
echo >> ${outd}/verkko.yml "ruk_type:            '${ruk_type}'"
echo >> ${outd}/verkko.yml "ruk_fract:           '${ruk_fract}'"
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
echo >> ${outd}/verkko.yml "#  This is the end."

#
#  If cnspaths (and cnsassembly) are defined and exist, set up to run just
#  consensus on the paths in ${cnspaths}, using various bits from an existing
#  Verkko assembly in ${cnsassembly}.
#

target="verkko"

if [ "x$cnspaths" != "x" ] ; then
    target="cnspath"

    cnspaths=$(fullpath $cnspaths)         #  Convert to full
    cnsassembly=$(fullpath $cnsassembly)   #  path.

    if [ ! -e $cnspaths ] ; then
        echo "Can't find --paths ${cnspaths}."
        exit
    fi

    if [ ! -d $cnsassembly ] ; then
        echo "Can't find --assembly ${cnsassembly}."
        exit
    fi

    #  Copy pieces from the previous assembly to the new run directory.  This
    #  is done - instead of symlinking - to prevent Snakemake from
    #  'accidentally' obliterating precious original files.

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

echo  > ${outd}/snakemake.sh "#!/bin/sh"
echo >> ${outd}/snakemake.sh ""
echo >> ${outd}/snakemake.sh "echo \"Launching $version ...\""
echo >> ${outd}/snakemake.sh ""
echo >> ${outd}/snakemake.sh "snakemake ${target} --nocolor \\"
echo >> ${outd}/snakemake.sh "  --directory . \\"
echo >> ${outd}/snakemake.sh "  --snakefile ${verkko}/Snakefile \\"
echo >> ${outd}/snakemake.sh "  --configfile verkko.yml \\"
echo >> ${outd}/snakemake.sh "  --reason \\"
echo >> ${outd}/snakemake.sh "  --keep-going \\"
echo >> ${outd}/snakemake.sh "  --rerun-incomplete \\"
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
