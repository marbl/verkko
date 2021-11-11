#!/bin/bash
set -eou

if [ "$#" -lt 2 ]; then
    #unitig-unrolled-ont-resolved.gfa unitig-popped-unitig-normal-connected-tip.gfa
    echo "Usage: $0 <in.gfa> <out.gfa>"
    exit 1
fi

SCRIPTS=$(dirname $(readlink -e $0))/../scripts

in_gfa=$1
out_gfa=$2

UntipRelative 30000 30000 0.1 0.1 < $in_gfa > connected-tip.gfa
$SCRIPTS/unitigify.py "utig3-" unitig-mapping-3.txt < connected-tip.gfa > unitig-normal-connected-tip.gfa
$SCRIPTS/pop_bubbles_keep_longest.py 10 < unitig-normal-connected-tip.gfa > popped-unitig-normal-connected-tip.gfa
$SCRIPTS/unitigify.py "utig4-" unitig-mapping-4.txt < popped-unitig-normal-connected-tip.gfa > $out_gfa
