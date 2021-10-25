#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    echo "launch.sh <config.yaml> <out_dir> [additional snakemake arguments]*"
    exit 239
fi

config=$(readlink -e $1)
mkdir -p $2
dir=$(readlink -e $2)

cd $(dirname $0)

export ROOT_PATH=$(readlink -e .)
export SCRIPT_PATH=$ROOT_PATH/../scripts
#Relies on having seqtk, winnowmap, MGB, GraphAigner (and UntipRelative) in PATH
module load seqtk
export PATH=/home/nurks2/git/Winnowmap/bin:/home/nurks2/git/MBG/bin:home/nurks2/git/GraphAligner/bin:$PATH

# -o /home/snurk/logs/ -e /home/snurk/logs/
echo "Will use config file: \"$config\""
echo "Output folder: \"$dir\""
echo "Additional snakemake arguments: \"${@:3}\""

snakemake --latency-wait 60 -j 4 --configfile $config --directory $dir "${@:3}"

cd -
