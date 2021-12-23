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

$SCRIPTS/untip_relative.py 35000 35000 0.1 0.1 < $in_gfa > connected-tip.gfa
$SCRIPTS/unitigify.py "utig3-" unitig-mapping-3.txt < connected-tip.gfa > unitig-normal-connected-tip.gfa
cat *mapping* > combined-nodemap-1.txt
rm -f combined-edges-1.gfa
cat hifi-resolved.gfa gapped-once-hifi-resolved.gfa gapped-twice-hifi-resolved.gfa gapped-hifi-resolved.gfa unrolled-hifi-resolved.gfa unitig-unrolled-hifi-resolved.gfa gapped-unitig-unrolled-hifi-resolved.gfa connected.gfa normal-connected.gfa ont-resolved-graph.gfa unrolled-ont-resolved.gfa unitig-unrolled-ont-resolved.gfa connected-tip.gfa unitig-normal-connected-tip.gfa | grep -P '^L' > combined-edges-1.gfa
grep -P '^S' hifi-resolved.gfa gapped-once-hifi-resolved.gfa gapped-twice-hifi-resolved.gfa gapped-hifi-resolved.gfa unrolled-hifi-resolved.gfa unitig-unrolled-hifi-resolved.gfa gapped-unitig-unrolled-hifi-resolved.gfa connected.gfa normal-connected.gfa ont-resolved-graph.gfa unrolled-ont-resolved.gfa unitig-unrolled-ont-resolved.gfa connected-tip.gfa unitig-normal-connected-tip.gfa | grep -v '\*' | awk '{print $2 "\t" length($3);}' > nodelens-1.txt
$SCRIPTS/get_original_coverage.py unitig-normal-connected-tip.gfa combined-nodemap-1.txt combined-edges-1.gfa nodelens-1.txt hifi_nodecov.csv > nodecov_hifi_fix.csv
$SCRIPTS/pop_bubbles_coverage_based.py nodecov_hifi_fix.csv < unitig-normal-connected-tip.gfa > popped-unitig-normal-connected-tip.gfa 2> popinfo.txt
$SCRIPTS/unitigify.py "utig4-" unitig-mapping-4.txt < popped-unitig-normal-connected-tip.gfa > $out_gfa
