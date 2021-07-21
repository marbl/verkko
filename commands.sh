MBG -i hifi.fa -o graph-multirle-nonblunt.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=msat --output-sequence-paths paths.gaf

scripts/calculate_coverage.py graph-multirle-nonblunt.gfa < paths.gaf > nodecovs-hifi.csv
scripts/insert_aln_gaps.py graph-multirle-nonblunt.gfa 2 < paths.gaf > gapped-graph-multirle-nonblunt.gfa
cut -f 6 < paths.gaf | scripts/unroll_tip_loops.py gapped-graph-multirle-nonblunt.gfa 3 > unrolled-graph-multirle-nonblunt.gfa

scripts/add_fake_alignments.py graph-multirle-nonblunt.gfa unrolled-graph-multirle-nonblunt.gfa paths.gaf nodecovs-hifi.csv fake-hifi-alns.gaf fake-hifi-nodecovs.csv 5
scripts/remove_small_tips.py unrolled-graph-multirle-nonblunt.gfa fake-hifi-nodecovs.csv fake-hifi-alns.gaf 2 7500 3 10 > untippedone-graph.gfa
scripts/remove_small_tips.py untippedone-graph.gfa fake-hifi-nodecovs.csv fake-hifi-alns.gaf 2 7500 5 20 > untippedtwo-graph.gfa
scripts/remove_lowcov_wrong_bubbles.py untippedtwo-graph.gfa fake-hifi-nodecovs.csv 3 7500 10 > untipped-graph.gfa
scripts/add_fake_alignments.py graph-multirle-nonblunt.gfa untipped-graph.gfa paths.gaf nodecovs-hifi.csv fake-hifi-alns.gaf fake-hifi-nodecovs.csv 5

scripts/estimate_unique_local.py untipped-graph.gfa paths.gaf 200000 10 > unique_nodes_hifi.txt
scripts/get_existing_paths.py untipped-graph.gfa < paths.gaf > paths.txt
scripts/find_bridges.py unique_nodes_hifi.txt < paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | scripts/remove_wrong_connections_2.py | sort > bridging_seq_all.txt
scripts/pick_majority_bridge.py < bridging_seq_all.txt > bridging_seq_picked_all.txt
cp bridging_seq_picked_all.txt bridging_seq_picked.txt
scripts/forbid_unbridged_tangles.py unique_nodes_hifi.txt untipped-graph.gfa bridging_seq_all.txt paths.txt nodecovs-hifi.csv 10 > forbidden_ends.txt
scripts/connect_uniques.py untipped-graph.gfa forbidden_ends.txt bridging_seq_picked.txt > hifi_connected.gfa
scripts/merge_unresolved_dbg_nodes.py < hifi_connected.gfa > normal-hifi_connected.gfa

scripts/add_fake_alignments.py graph-multirle-nonblunt.gfa normal-hifi_connected.gfa paths.gaf nodecovs-hifi.csv fake-hifi-alns.gaf fake-hifi-nodecovs.csv 5
/usr/bin/time -v scripts/resolve_triplets_kmerify.py normal-hifi_connected.gfa fake-hifi-paths.txt fake-hifi-nodecovs.csv 15000 2 2 < fake-hifi-alns.gaf > hifi-resolved-graph.gfa 2> stderr_hifi_resolved_graph.txt

grep -P '^S\t' < hifi-resolved-graph.gfa | awk 'BEGIN{print "node\tlength\tcoverage";}{if (length($3) > 15000) {print $2 "\t" length($3) "\t" 10;} else {print $2 "\t" length($3) "\t" 0;}}' > fake-coverages.csv
scripts/translate_uniques.py hifi-resolved-graph.gfa < unique_nodes_hifi.txt > unique_nodes_hifi_translate.txt
scripts/find_bridges.py unique_nodes_hifi_translate.txt < fake-hifi-paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | scripts/remove_wrong_connections_2.py | sort > bridging_seq_all.txt
scripts/pick_majority_bridge.py < bridging_seq_all.txt > bridging_seq_picked_all.txt
cp bridging_seq_picked_all.txt bridging_seq_picked.txt
scripts/forbid_unbridged_tangles.py unique_nodes_hifi_translate.txt hifi-resolved-graph.gfa bridging_seq_all.txt fake-hifi-paths.txt fake-coverages.csv 10 > forbidden_ends.txt
scripts/connect_uniques.py hifi-resolved-graph.gfa forbidden_ends.txt bridging_seq_picked.txt > hifi_connected_twice.gfa
scripts/merge_unresolved_dbg_nodes.py < hifi_connected_twice.gfa > normal-hifi_connected_twice.gfa
scripts/unitigify.py < normal-hifi_connected_twice.gfa > unitig-normal-hifi_connected_twice.gfa

scripts/unroll_tip_loops.py unitig-normal-hifi_connected_twice.gfa 3 < fake-hifi-paths.txt > unrolled-hifi-resolved.gfa
scripts/unitigify.py < unrolled-hifi-resolved.gfa > unitig-unrolled-hifi-resolved.gfa

UntipRelative 15000 15000 0.1 0.1 < unitig-unrolled-hifi-resolved.gfa > hifi-resolved-graph-tip.gfa
scripts/pop_bubbles_keep_longest.py 10 < hifi-resolved-graph-tip.gfa > popped-hifi-resolved-graph-tip.gfa
scripts/unitigify.py < popped-hifi-resolved-graph-tip.gfa > unitig-popped-hifi-resolved-graph-tip.gfa

GraphAligner -t 32 -g unitig-unrolled-hifi-resolved.gfa -f ont.fa -a alns-ont.gaf --seeds-mxm-length 30 --seeds-mem-count 10000 -b 15 --multimap-score-fraction 0.99 --precise-clipping 0.7 --min-alignment-score 5000 --hpc-collapse-reads --discard-cigar 1> stdout_ga_ont.txt 2> stderr_ga_ont.txt

awk -F '\t' '{if ($4-$3 >= $2*0.8 && $12 >= 20) print;}' < alns-ont.gaf > alns-ont-filter.gaf
scripts/trim_dbg_alignment.py unitig-unrolled-hifi-resolved.gfa 1500 < alns-ont-filter.gaf > alns-ont-filter-trim.gaf
scripts/calculate_coverage.py unitig-unrolled-hifi-resolved.gfa < alns-ont-filter-trim.gaf > nodecovs-ont.csv

cut -f 6 < alns-ont-filter-trim.gaf > paths.txt
awk -F '\t' '{if ($12 >= 20) print;}' < alns-ont.gaf > alns-ont-mapqfilter.gaf
scripts/insert_aln_gaps.py unitig-unrolled-hifi-resolved.gfa 3 < alns-ont-mapqfilter.gaf > gapped-unitig-unrolled-hifi-resolved.gfa
awk '{if ($2 >= 150000) {sum += $2*$3; count += $2;}}END{print sum/count;}' < nodecovs-ont.csv
scripts/estimate_unique_local.py gapped-unitig-unrolled-hifi-resolved.gfa alns-ont-filter-trim.gaf 200000 30 > unique_nodes_ont_coverage.txt
# scripts/translate_uniques.py normal-hifi_connected_twice.gfa < unique_nodes_hifi.txt > translated_uniques.txt
# scripts/translate_nodes_by_seq.py normal-hifi_connected_twice.gfa unitig-unrolled-hifi-resolved.gfa < translated_uniques.txt > unique_nodes_ont_translated.txt
# cat unique_nodes_ont_coverage.txt unique_nodes_ont_translated.txt | sort | uniq > unique_nodes_ont.txt
cp unique_nodes_ont_coverage.txt unique_nodes_ont.txt

scripts/find_bridges.py unique_nodes_ont.txt < paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | scripts/remove_wrong_connections_2.py | sort > bridging_seq_all.txt
scripts/pick_majority_bridge.py < bridging_seq_all.txt > bridging_seq_picked_all.txt
cp bridging_seq_picked_all.txt bridging_seq_picked.txt
scripts/forbid_unbridged_tangles.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridging_seq_all.txt paths.txt nodecovs-ont.csv 30 > forbidden_ends.txt
scripts/connect_uniques.py gapped-unitig-unrolled-hifi-resolved.gfa forbidden_ends.txt bridging_seq_picked.txt > connected.gfa

scripts/merge_unresolved_dbg_nodes.py < connected.gfa > normal-connected.gfa
scripts/add_fake_alignments.py unitig-unrolled-hifi-resolved.gfa normal-connected.gfa alns-ont-filter-trim.gaf nodecovs-ont.csv fake-ont-alns.gaf fake-ont-nodecovs.csv 10
/usr/bin/time -v scripts/resolve_triplets_kmerify.py normal-connected.gfa fake-ont-paths.txt fake-ont-nodecovs.csv 100000 3 5 3 2 < fake-ont-alns.gaf > ont-resolved-graph.gfa 2> stderr_ont_resolved_graph.txt

scripts/unroll_tip_loops.py ont-resolved-graph.gfa 3 < fake-ont-paths.txt > unrolled-ont-resolved.gfa
scripts/unitigify.py < unrolled-ont-resolved.gfa > unitig-unrolled-ont-resolved.gfa



UntipRelative 30000 30000 0.1 0.1 < ont-resolved-graph.gfa > connected-tip.gfa
scripts/unitigify.py < connected-tip.gfa > unitig-normal-connected-tip.gfa
scripts/pop_bubbles_keep_longest.py 10 < unitig-normal-connected-tip.gfa > popped-unitig-normal-connected-tip.gfa
scripts/unitigify.py < popped-unitig-normal-connected-tip.gfa > unitig-popped-unitig-normal-connected-tip.gfa
