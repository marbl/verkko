#!/bin/sh
#
#  Verkko, Finnish.  Net, network, mesh, web, grid, grill, fishnet, network, graph.
#
#  Requires:
#    seqtk
#    winnowmap
#    MGB
#    GraphAligner and UntipRelative
#    seqrequester
#
#  profiles - searches for folder in /etc/xdg/snakemake and ~/.config/snakemake,
#  or absolute or relative path.  Folder must contain config.yaml.
#    --cluster qsub ==> cluster: qsub
#  or set via env var SNAKEMAKE_PROfiLE
#

wd=`pwd`

hifi=""
nano=""
outd=""

errors=""
verkko=$VERKKO

if [ "x$verkko" = "x" ] ; then   #  Set verkko to a bogus value
  verkko='$VERKKO'               #  for error reporting.
fi

mbg=
graphaligner=
python=

snakeopts="--unlock"
snakeopts="--cleanup-metadata 3-align/split"
snakeopts="--cleanup-metadata 6-layout/layout.txt --cleanup-metadata 6-layout/gaps.txt"
snakeopts=""

#  Algorithm parameters.

#  buildGraph, parameters for MBG
mbg_baseK=1001
mbg_maxK=15000
mbg_window=100

#  split_ont, partitioning ont reads for alignment
spl_bases=3000000000
spl_reads=150000
spl_min_length=0

#  align_ont, alignment of ont reads to the initial graph
ali_mxm_length=30
ali_mem_count=1000
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

#  Run parameters.

#  build-graph
mbg_n_cpus=4
mbg_mem_gb=128
mbg_time_h=24

#  process_graph
utg_n_cpus=1
utg_mem_gb=64
utg_time_h=24

#  split_ont
spl_n_cpus=1
spl_mem_gb=8
spl_time_h=24

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

if [ $# -eq 0 ] ; then
  echo "usage: $0 -d <work-directory> [options] -hifi reads.fastq.gz ... -nano reads.fastq.gz ..."
  exit
fi


while [ $# -gt 0 ] ; do
    opt=$1
    arg=$2

    shift

    #
    #  Run options
    #

    if   [ "$opt" = "-d" ] ; then
        outd=$arg
        shift

    elif [ "$opt" = "--python" ] ; then
        python=$arg
        shift

    elif [ "$opt" = "--mbg" ] ; then
        mbg=$arg
        shift

    elif [ "$opt" = "--graphaligner" ] ; then
        graphaligner=$arg
        shift

    elif [ "$opt" = "--profile" ] ; then
        profile=$arg
        shift

    elif [ "$opt" = "--snakeopts" ] ; then
        snakeopts=$arg
        shift

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
    #  MBG options
    #

    elif [ "$opt" = "--base-k" ] ; then
        mbg_baseK=$arg
        shift

    elif [ "$opt" = "--max-k" ] ; then
        mbg_maxK=$arg
        shift

    elif [ "$opt" = "--window" ] ; then
        mbg_window=$arg
        shift

    elif [ "$opt" = "--threads" ] ; then
        mbg_threads=$arg
        shift

    #
    #  splitONT options
    #

    elif [ "$opt" = "--split-bases" ] ; then
        spl_bases=$arg
        shift

    elif [ "$opt" = "--split-reads" ] ; then
        spl_reads=$arg
        shift

    elif [ "$opt" = "--min-ont-length" ] ; then
        spl_min_length=$arg
        shift

    #
    #  alignONT options
    #

    elif [ "$opt" = "--seed-min-length" ] ; then
        ali_mxm_length=$arg
        shift

    elif [ "$opt" = "--seed-max-length" ] ; then
        ali_mem_count=$arg
        shift

    elif [ "$opt" = "--align-bandwidth" ] ; then
        ali_bandwidth=$arg
        shift

    elif [ "$opt" = "--score-fraction" ] ; then
        ali_multi_score_f=$arg
        shift

    elif [ "$opt" = "--min-identity" ] ; then
        ali_clipping=$arg
        shift

    elif [ "$opt" = "--min-score" ] ; then
        ali_min_score=$arg
        shift

    elif [ "$opt" = "--end-clipping" ] ; then
        ali_end_clipping=$arg
        shift

    elif [ "$opt" = "--incompatible-cutoff" ] ; then
        ali_incompat_cutoff=$arg
        shift

    elif [ "$opt" = "--max-traces" ] ; then
        ali_max_trace=$arg
        shift

    else
        errors="${errors}Unknown option '$opt'.\n"
    fi
done

#
#  Set stuff not set by options.
#

if [ "x$mbg" = "x" ] ; then
  mbg=${verkko}/bin/MBG
fi

if [ "x$graphaligner" = "x" ] ; then
  graphaligner=${verkko}/bin/GraphAligner
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

if [ ! -e "$verkko/verkko.sh" ] ; then
    errors="${errors}Can't find 'verkko.sh' in directory VERKKO = '$verkko/'.\n"
fi

if [ ! -e "$verkko/scripts/get_layout_from_mbg.py" ] ; then
    errors="${errors}Can't find Verkko scripts/ directory at '$verkko/scripts/'.\n"
fi

if [ ! -e "$mbg" ] ; then
    errors="${errors}Can't find MBG executable at '$mbg'.\n"
fi

if [ ! -e "$graphaligner" ] ; then
    errors="${errors}Can't find GraphAligner executable at '$graphaligner'.\n"
fi


if [ "x$errors" != "x" ] ; then
    echo "ERRORS."
    printf "${errors}"
    exit 0
fi


mkdir -p ${outd}
cd       ${outd}

#echo "Will use config file:           '$conf'"
echo "Output folder:                  '$outd'"
echo "Additional snakemake arguments: '$snakeopts'"

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
#echo >> verkko.yml "HIFI_READS:          '[ ${hifi} ]'"
#echo >> verkko.yml "ONT_READS:           '[ ${nano} ]'"

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
echo >> verkko.yml "#  build-graph, MBG"
echo >> verkko.yml "mbg_baseK:           '${mbg_baseK}'"
echo >> verkko.yml "mbg_maxK:            '${mbg_maxK}'"
echo >> verkko.yml "mbg_window:          '${mbg_window}'"
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

export PATH=${verkko}/bin:${verkko}/scripts:$PATH

#  --cluster-config cluster.json \
#  --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" \

#  jobs  - N cloud/cluster jobs in parallel (for local, an alais for -cores)
#  cores - N total number of cores used over all jobs

#echo >> snakemake.sh "  --jobs 60 \\"
#echo >> snakemake.sh "  --verbose \\"
#echo >> snakemake.sh "  --cleanup-metadata 6-layout/layout.txt \\"

#  For interacrtive
#echo >> snakemake.sh "  --cores 4 \\"
#echo >> snakemake.sh "  --profile ${verkko}/profiles \\"

echo  > snakemake.sh "#!/bin/sh"
echo >> snakemake.sh ""
echo >> snakemake.sh "snakemake --nocolor \\"
echo >> snakemake.sh "  --directory . \\"
echo >> snakemake.sh "  --snakefile ${verkko}/Snakefile \\"
echo >> snakemake.sh "  --latency-wait 60 \\"
echo >> snakemake.sh "  --jobs 1000 \\"
echo >> snakemake.sh "  --cores 16 \\"
echo >> snakemake.sh "  --configfile verkko.yml \\"
#echo >> snakemake.sh "  --printshellcmds \\"
echo >> snakemake.sh "  --reason \\"
echo >> snakemake.sh "  ${snakeopts}"
echo >> snakemake.sh ""

chmod +x snakemake.sh

./snakemake.sh

exit 0
