#!/bin/sh
set -eou

if [ "$#" -lt 5 ]; then
    #unitig-unrolled-hifi-resolved.gfa alns-ont.gaf unitig-unrolled-ont-resolved.gfa 5 20 10 5
    echo "Usage: $0 <in.gfa> <ont_paths.gaf> <out.gfa> <min_allowed_coverage> <resolve_steps>+"
    exit 1
fi

SCRIPTS=$(dirname $(readlink -e $0))/../scripts

in_gfa=$1
ont_paths=$2
out_gfa=$3
#for triplet resolution
min_allowed_coverage=$4
echo "In triplet resolution will use min_allowed_coverage of $min_allowed_coverage and resolve_steps ${@:5}"

awk -F '\t' '{if ($4-$3 >= $2*0.8 && $12 >= 20) print;}' < $ont_paths > alns-ont-filter.gaf
$SCRIPTS/trim_dbg_alignment.py $in_gfa 1500 < alns-ont-filter.gaf > alns-ont-filter-trim.gaf
$SCRIPTS/calculate_coverage.py $in_gfa < alns-ont-filter-trim.gaf > nodecovs-ont.csv

cut -f 6 < alns-ont-filter-trim.gaf > paths.txt
awk -F '\t' '{if ($12 >= 20) print;}' < $ont_paths > alns-ont-mapqfilter.gaf
$SCRIPTS/insert_aln_gaps.py $in_gfa 3 50 gaps-ont.gaf gapont < alns-ont-mapqfilter.gaf > gapped-unitig-unrolled-hifi-resolved.gfa
awk '{if ($2 >= 100000) {sum += $2*$3; count += $2;}}END{print sum/count;}' < nodecovs-ont.csv

cat *mapping* > combined-nodemap-uniques.txt
rm -f combined-edges-uniques.gfa
cat *.gfa | grep -P '^L' > combined-edges-uniques.gfa
grep -P '^S' *.gfa | grep -v '\*' | awk '{print $2 "\t" length($3);}' > nodelens-uniques.txt
$SCRIPTS/get_original_coverage.py gapped-unitig-unrolled-hifi-resolved.gfa combined-nodemap-uniques.txt combined-edges-uniques.gfa nodelens-uniques.txt hifi_nodecov.csv > hifi-nodecov-gapped-unitig-unrolled-hifi-resolved.csv

$SCRIPTS/estimate_unique_local.py gapped-unitig-unrolled-hifi-resolved.gfa hifi-nodecov-gapped-unitig-unrolled-hifi-resolved.csv alns-ont-filter-trim.gaf 100000 30 0.8 > unique_nodes_ont_coverage.txt
# $SCRIPTS/translate_uniques.py normal-hifi_connected_twice.gfa < unique_nodes_hifi.txt > translated_uniques.txt
# $SCRIPTS/translate_nodes_by_seq.py normal-hifi_connected_twice.gfa $in_gfa < translated_uniques.txt > unique_nodes_ont_translated.txt
# cat unique_nodes_ont_coverage.txt unique_nodes_ont_translated.txt | sort | uniq > unique_nodes_ont.txt
$SCRIPTS/fix_diploid_unique_nodes.py unique_nodes_ont_coverage.txt nodecovs-ont.csv gapped-unitig-unrolled-hifi-resolved.gfa > unique_nodes_diploidfix.txt
cp unique_nodes_diploidfix.txt unique_nodes_ont.txt

$SCRIPTS/find_bridges.py unique_nodes_ont.txt < paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | $SCRIPTS/remove_wrong_connections_2.py forbidden_wrong_connections.txt | sort > bridging_seq_all.txt
$SCRIPTS/pick_majority_bridge.py forbidden_minority_bridges.txt < bridging_seq_all.txt > bridging_seq_picked_all.txt
$SCRIPTS/remove_crosslink_paths.py unique_nodes_ont.txt bridging_seq_picked_all.txt bridges.txt > bridges_fixcrosslink.txt 2> forbidden_crosslinks.txt
$SCRIPTS/fix_diploid_paths.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridges_fixcrosslink.txt bridges.txt 3 > bridging_seq_diploidfix_all.txt
cp bridging_seq_diploidfix_all.txt bridging_seq_picked.txt
# forbidden_wrong_connections.txt deliberately not included here so that if that causes a gap, the tangle is forbidden
cat forbidden_crosslinks.txt forbidden_minority_bridges.txt > bridging_seq_forbidden.txt
$SCRIPTS/forbid_unbridged_tangles.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridging_seq_forbidden.txt bridging_seq_picked.txt paths.txt nodecovs-ont.csv 30 > forbidden_ends.txt
$SCRIPTS/connect_uniques.py gapped-unitig-unrolled-hifi-resolved.gfa forbidden_ends.txt bridging_seq_picked.txt > connected.gfa

$SCRIPTS/merge_unresolved_dbg_nodes.py < connected.gfa > normal-connected.gfa
$SCRIPTS/get_bridge_mapping.py normal-connected.gfa gapped-unitig-unrolled-hifi-resolved.gfa > bridge_mapping.txt
$SCRIPTS/add_fake_alignments.py $in_gfa normal-connected.gfa alns-ont-filter-trim.gaf nodecovs-ont.csv fake-ont-alns.gaf fake-ont-nodecovs.csv 10
$SCRIPTS/add_fake_bridging_paths.py forbidden_ends.txt bridging_seq_picked.txt fake-ont-nodecovs-once.csv fake-ont-nodecovs.csv 10 >> fake-ont-alns.gaf
#FIXME parameterize
/usr/bin/time -v $SCRIPTS/resolve_triplets_kmerify.py normal-connected.gfa fake-ont-paths.txt fake-ont-nodecovs.csv resolve-mapping.txt 100000 $min_allowed_coverage ${@:5} < fake-ont-alns.gaf > ont-resolved-graph.gfa 2> stderr_ont_resolved_graph.txt

$SCRIPTS/unroll_tip_loops.py ont-resolved-graph.gfa 3 < fake-ont-paths.txt > unrolled-ont-resolved.gfa
$SCRIPTS/get_unroll_mapping.py ont-resolved-graph.gfa unrolled-ont-resolved.gfa > unroll_mapping_2.txt
$SCRIPTS/unitigify.py "utig2-" unitig-mapping-2.txt < unrolled-ont-resolved.gfa > $out_gfa
