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

#  Verkko, Finnish.  Net, network, mesh, web, grid, grill, fishnet, network, graph.
#

version=""
help=""

hifi=""
nano=""
outd=""

errors=""
verkko=""

mbg=""
graphaligner=""
python=

grid="local"
local_cpus=all
local_mem=64

snakeopts=""


#  If environment variable VERKKO is set, assume that is the path to our
#  installed directory, otherwise, figure it out from the shell script name
#  and set VERKKO.
#
#  This is to allow submitting verkko.sh directly to the grid; the grid will
#  (typically) copy the shell script to a spool directory and run it there -
#  so any auto-detected path based on the script name will be incorrect.

if [ -e "${VERKKO}/verkko.sh" ] ; then
  verkko=$VERKKO
else
  echo $0 | grep -q ^/            #  If a relative path, prepend
  if [ $? -ne 0 ] ; then          #  PWD to make it an absolute path.
    f=$0

    echo $0 | grep -q ^\\./       #  If the path looks like ./ something
    if [ $? -eq 0 ] ;then         #  remove the start ./ from the command
       f=`echo $0 | cut -d'/' -f2-`
    fi
    verkko=`dirname $PWD/$f`
  else
    verkko=`dirname $0`
  fi

  #  If we're running out of the src/ directory, we need to make
  #  a symlink to the Canu binaries, but otherwise, the Verkko
  #  components are in the correct place.

  if   [ -e "${verkko}/verkko.sh" ] ; then
    if [ ! -e "bin" ] ; then
      cd ${verkko} && ln -s build/bin .
    fi

  #  Otherwise, massage the path to point to lib/verkko instead of bin/,

  else
    verkko=`dirname $verkko`
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
mbg_max_resolution=2000

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
ali_mem_gb=160
ali_time_h=24

#  process_ont_paths
pop_n_cpus=1
pop_mem_gb=64
pop_time_h=24

#  untip
utp_n_cpus=1
utp_mem_gb=64
utp_time_h=24

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
par_mem_gb=120
par_time_h=24

#  cns
cns_n_cpus=8
cns_mem_gb=200
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

    elif [ "$opt" = "-d" ] ;             then outd=$arg;          shift
    elif [ "$opt" = "--python" ] ;       then python=$arg;        shift
    elif [ "$opt" = "--mbg" ] ;          then mbg=$arg;           shift
    elif [ "$opt" = "--graphaligner" ] ; then graphaligner=$arg;  shift
    elif [ "$opt" = "--local" ] ;        then grid="local";
    elif [ "$opt" = "--sge" ] ;          then grid="slurm-sge";
    elif [ "$opt" = "--slurm" ] ;        then grid="slurm-sge";
    elif [ "$opt" = "--lsf" ] ;          then grid="lsf";
    elif [ "$opt" = "--local-memory" ] ; then local_mem=$arg;     shift
    elif [ "$opt" = "--local-cpus" ] ;   then local_cpus=$arg;    shift
    elif [ "$opt" = "--snakeopts" ] ;    then snakeopts=$arg;     shift

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
            #errors="${errors}Can\'t find HiFi reads file \'$arg\'.\n"
        done

    elif [ "$opt" = "--nano" ] ; then
        while [ -e "$arg" ] ; do
            if [ -e "/$arg" ] ; then
                nano="$nano $arg"
            else
                nano="$nano `pwd`/$arg"
            fi
            shift
            arg=$1
            #errors="${errors}Can\'t find Nanopore reads file \'$arg\'.\n"
        done

    #
    #
    #

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
    elif [ "$opt" = "--red-run" ] ;  then red_n_cpus=$1; red_mem_gb=$2; red_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--mbg-run" ] ;  then mbg_n_cpus=$1; mbg_mem_gb=$2; mbg_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--utg-run" ] ;  then utg_n_cpus=$1; utg_mem_gb=$2; utg_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--spl-run" ] ;  then spl_n_cpus=$1; spl_mem_gb=$2; spl_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--ali-run" ] ;  then ali_n_cpus=$1; ali_mem_gb=$2; ali_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--pop-run" ] ;  then pop_n_cpus=$1; pop_mem_gb=$2; pop_time_h=$3; shift; shift; shift;
    elif [ "$opt" = "--utp-run" ] ;  then utp_n_cpus=$1; utp_mem_gb=$2; utp_time_h=$3; shift; shift; shift;
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

if [ "x$graphaligner" = "x" ] ; then
  graphaligner=${verkko}/bin/GraphAligner
fi
if [ ! -e $graphaligner ] ; then
  graphaligner=$(which GraphAligner)
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
    errors="${errors}No PacBio HiFi reads (-hifi) supplied.\n"
fi

if [ "x$nano" = "x" ] ; then
    errors="${errors}No Oxford Nanopore reads (-nano) supplied.\n"
fi

#           bin/seqrequester

for exe in bin/findErrors \
           bin/fixErrors \
           bin/layoutToPackage \
           bin/meryl \
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
    echo "usage: $0 -d <output-directory> -hifi <hifi-reads ...> -nano <nanopore-reads ...>"
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
    echo "    --no-correction"
    echo ""
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
#    Make a work directory for us.
#    Write a yaml config file for snakemake.
#    Create a script to run snakemake.
#

mkdir -p ${outd}
cd       ${outd}

echo  > verkko.yml "#  Generated automatically by verkko.sh."
echo >> verkko.yml "#  Changes will be overwritten."
echo >> verkko.yml ""
echo >> verkko.yml "VERKKO:              '${verkko}'"
echo >> verkko.yml ""
echo >> verkko.yml "MBG:                 '${mbg}'"
echo >> verkko.yml "GRAPHALIGNER:        '${graphaligner}'"
echo >> verkko.yml ""
echo >> verkko.yml "PYTHON:              '${python}'"
echo >> verkko.yml ""
echo >> verkko.yml "HIFI_READS:"
for h in ${hifi} ; do
  echo >> verkko.yml " - '$h'"
done
echo >> verkko.yml "ONT_READS:"
for o in ${nano} ; do
  echo >> verkko.yml " - '$o'"
done
echo >> verkko.yml ""
echo >> verkko.yml "#  Algorithm parameters."
echo >> verkko.yml ""
echo >> verkko.yml "#  buildStore, countKmers and computeOverlaps"
echo >> verkko.yml "correction_enabled:  '${correction_enabled}'"
echo >> verkko.yml "mer_size:            '${mer_size}'"
echo >> verkko.yml "mer_threshold:       '${mer_threshold}'"
echo >> verkko.yml ""
echo >> verkko.yml "cor_min_read:        '${cor_min_read}'"
echo >> verkko.yml "cor_min_overlap:     '${cor_min_overlap}'"
echo >> verkko.yml "cor_hash_bits:       '${cor_hash_bits}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  build-graph, MBG"
echo >> verkko.yml "mbg_baseK:           '${mbg_baseK}'"
echo >> verkko.yml "mbg_maxK:            '${mbg_maxK}'"
echo >> verkko.yml "mbg_window:          '${mbg_window}'"
echo >> verkko.yml "mbg_max_resolution:  '${mbg_max_resolution}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  split_ont"
echo >> verkko.yml "spl_bases:           '${spl_bases}'"
echo >> verkko.yml "spl_reads:           '${spl_reads}'"
echo >> verkko.yml "spl_min_length:      '${spl_min_length}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  align_ont, GraphAligner"
echo >> verkko.yml "ali_mxm_length:      '${ali_mxm_length}'"
echo >> verkko.yml "ali_mem_count:       '${ali_mem_count}'"
echo >> verkko.yml "ali_bandwidth:       '${ali_bandwidth}'"
echo >> verkko.yml "ali_multi_score_f:   '${ali_multi_score_f}'"
echo >> verkko.yml "ali_clipping:        '${ali_clipping}'"
echo >> verkko.yml "ali_min_score:       '${ali_min_score}'"
echo >> verkko.yml "ali_end_clipping:    '${ali_end_clipping}'"
echo >> verkko.yml "ali_incompat_cutoff: '${ali_incompat_cutoff}'"
echo >> verkko.yml "ali_max_trace:       '${ali_max_trace}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  process_ont_paths"
echo >> verkko.yml "pop_min_allowed_cov: '${pop_min_allowed_cov}'"
echo >> verkko.yml "pop_resolve_steps:   '${pop_resolve_steps}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  Run parameters."
echo >> verkko.yml ""
echo >> verkko.yml "#  buildStore, countKmers and computeOverlaps"
echo >> verkko.yml "sto_n_cpus:          '${sto_n_cpus}'"
echo >> verkko.yml "sto_mem_gb:          '${sto_mem_gb}'"
echo >> verkko.yml "sto_time_h:          '${sto_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "mer_n_cpus:          '${mer_n_cpus}'"
echo >> verkko.yml "mer_mem_gb:          '${mer_mem_gb}'"
echo >> verkko.yml "mer_time_h:          '${mer_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "ovb_n_cpus:          '${ovb_n_cpus}'"
echo >> verkko.yml "ovb_mem_gb:          '${ovb_mem_gb}'"
echo >> verkko.yml "ovb_time_h:          '${ovb_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "red_n_cpus:          '${red_n_cpus}'"
echo >> verkko.yml "red_mem_gb:          '${red_mem_gb}'"
echo >> verkko.yml "red_time_h:          '${red_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  build-graph"
echo >> verkko.yml "mbg_n_cpus:          '${mbg_n_cpus}'"
echo >> verkko.yml "mbg_mem_gb:          '${mbg_mem_gb}'"
echo >> verkko.yml "mbg_time_h:          '${mbg_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  process_graph"
echo >> verkko.yml "utg_n_cpus:          '${utg_n_cpus}'"
echo >> verkko.yml "utg_mem_gb:          '${utg_mem_gb}'"
echo >> verkko.yml "utg_time_h:          '${utg_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  split_ont"
echo >> verkko.yml "spl_n_cpus:          '${spl_n_cpus}'"
echo >> verkko.yml "spl_mem_gb:          '${spl_mem_gb}'"
echo >> verkko.yml "spl_time_h:          '${spl_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  align_ont"
echo >> verkko.yml "ali_n_cpus:          '${ali_n_cpus}'"
echo >> verkko.yml "ali_mem_gb:          '${ali_mem_gb}'"
echo >> verkko.yml "ali_time_h:          '${ali_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  process_ont_paths"
echo >> verkko.yml "pop_n_cpus:          '${pop_n_cpus}'"
echo >> verkko.yml "pop_mem_gb:          '${pop_mem_gb}'"
echo >> verkko.yml "pop_time_h:          '${pop_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  untip"
echo >> verkko.yml "utp_n_cpus:          '${utp_n_cpus}'"
echo >> verkko.yml "utp_mem_gb:          '${utp_mem_gb}'"
echo >> verkko.yml "utp_time_h:          '${utp_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  create_layout"
echo >> verkko.yml "lay_n_cpus:          '${lay_n_cpus}'"
echo >> verkko.yml "lay_mem_gb:          '${lay_mem_gb}'"
echo >> verkko.yml "lay_time_h:          '${lay_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  get_ont_subset"
echo >> verkko.yml "sub_n_cpus:          '${sub_n_cpus}'"
echo >> verkko.yml "sub_mem_gb:          '${sub_mem_gb}'"
echo >> verkko.yml "sub_time_h:          '${sub_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  partition_consensus"
echo >> verkko.yml "par_n_cpus:          '${par_n_cpus}'"
echo >> verkko.yml "par_mem_gb:          '${par_mem_gb}'"
echo >> verkko.yml "par_time_h:          '${par_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  cns"
echo >> verkko.yml "cns_n_cpus:          '${cns_n_cpus}'"
echo >> verkko.yml "cns_mem_gb:          '${cns_mem_gb}'"
echo >> verkko.yml "cns_time_h:          '${cns_time_h}'"
echo >> verkko.yml ""
echo >> verkko.yml "#  This is the end."

echo  > snakemake.sh "#!/bin/sh"
echo >> snakemake.sh ""
echo >> snakemake.sh "echo \"Launching $version ...\""
echo >> snakemake.sh ""
echo >> snakemake.sh "snakemake --nocolor \\"
echo >> snakemake.sh "  --directory . \\"
echo >> snakemake.sh "  --snakefile ${verkko}/Snakefile \\"
echo >> snakemake.sh "  --configfile verkko.yml \\"
echo >> snakemake.sh "  --reason \\"
echo >> snakemake.sh "  --keep-going \\"
echo >> snakemake.sh "  --rerun-incomplete \\"
if [ $grid = "local" ] ; then
  echo >> snakemake.sh "  --latency-wait 2 \\"
  echo >> snakemake.sh "  --cores ${local_cpus} \\"
  echo >> snakemake.sh "  --resources mem_gb=${local_mem} \\"
else
  echo >> snakemake.sh "  --latency-wait 30 \\"
  echo >> snakemake.sh "  --jobs 1000 \\"
  echo >> snakemake.sh "  --profile ${verkko}/profiles \\"
  echo >> snakemake.sh "  --restart-times 1 \\"
  echo >> snakemake.sh "  --max-jobs-per-second 10 \\"
  echo >> snakemake.sh "  --max-status-checks-per-second 0.02 \\"
  echo >> snakemake.sh "  --local-cores 1 \\"
fi
echo >> snakemake.sh "  ${snakeopts}"
echo >> snakemake.sh ""

chmod +x snakemake.sh

#export PATH=${verkko}/bin:${verkko}/scripts:$PATH

./snakemake.sh

exit 0
