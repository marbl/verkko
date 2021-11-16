#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "launch.sh <config.yaml> <out_dir> [additional snakemake arguments]*"
    exit 239
fi

config=$(readlink -e $1)
mkdir -p $2
dir=$(readlink -e $2)

export LAUNCH_PATH=$(readlink -e $PWD)

cd $(dirname $0)

export ROOT_PATH=$(readlink -e .)
export SCRIPT_PATH=$ROOT_PATH/../scripts

#Relies on having seqtk, winnowmap, MGB, GraphAigner (and UntipRelative) in PATH
module load seqtk
export PATH=/home/nurks2/git/Winnowmap/bin:/data/rautiainenma/MBG/bin:/home/nurks2/git/GraphAligner/bin:/data/korens/devel/canu-verkko/build/bin:$PATH

# -o /home/snurk/logs/ -e /home/snurk/logs/
echo "Will use config file: \"$config\""
echo "Output folder: \"$dir\""
echo "Additional snakemake arguments: \"${@:3}\""

snakemake --latency-wait 60 -j 60 --local-cores 4 --cluster-config cluster.json --configfile $config --cluster "sbatch --ntasks 1 --mem {cluster.mem}G --cpus-per-task {cluster.n} --time {cluster.time}" --directory $dir "${@:3}"

cd -
