MBG -i hifi.fa -o hifi-resolved.gfa -k 1001 -w 100 -a 1 -u 2 --error-masking=collapse-msat --output-sequence-paths paths.gaf -r 15000

scripts/insert_aln_gaps.py hifi-resolved.gfa 2 50 /dev/null < paths.gaf > gapped-once-hifi-resolved.gfa
scripts/insert_aln_gaps.py gapped-once-hifi-resolved.gfa 2 300 /dev/null < paths.gaf > gapped-twice-hifi-resolved.gfa
scripts/insert_aln_gaps.py gapped-twice-hifi-resolved.gfa 1 5 /dev/null < paths.gaf > gapped-hifi-resolved.gfa
cut -f 6 < paths.gaf | scripts/unroll_tip_loops.py gapped-hifi-resolved.gfa 5 > unrolled-hifi-resolved.gfa
scripts/unitigify.py < unrolled-hifi-resolved.gfa > unitig-unrolled-hifi-resolved.gfa

# only for evaluating hifi-only CHM13 haploid assemblies
# UntipRelative 15000 15000 0.1 0.1 < unitig-unrolled-hifi-resolved.gfa > hifi-resolved-graph-tip.gfa
# scripts/pop_bubbles_keep_longest.py 10 < hifi-resolved-graph-tip.gfa > popped-hifi-resolved-graph-tip.gfa
# scripts/unitigify.py < popped-hifi-resolved-graph-tip.gfa > unitig-popped-hifi-resolved-graph-tip.gfa

GraphAligner -t 32 -g unitig-unrolled-hifi-resolved.gfa -f ont.fa -a alns-ont.gaf --seeds-mxm-length 30 --seeds-mem-count 10000 -b 15 --multimap-score-fraction 0.99 --precise-clipping 0.7 --min-alignment-score 5000 --hpc-collapse-reads --discard-cigar 1> stdout_ga_ont.txt 2> stderr_ga_ont.txt

awk -F '\t' '{if ($4-$3 >= $2*0.8 && $12 >= 20) print;}' < alns-ont.gaf > alns-ont-filter.gaf
scripts/trim_dbg_alignment.py unitig-unrolled-hifi-resolved.gfa 1500 < alns-ont-filter.gaf > alns-ont-filter-trim.gaf
scripts/calculate_coverage.py unitig-unrolled-hifi-resolved.gfa < alns-ont-filter-trim.gaf > nodecovs-ont.csv

cut -f 6 < alns-ont-filter-trim.gaf > paths.txt
awk -F '\t' '{if ($12 >= 20) print;}' < alns-ont.gaf > alns-ont-mapqfilter.gaf
scripts/insert_aln_gaps.py unitig-unrolled-hifi-resolved.gfa 3 50 used_ont.txt < alns-ont-mapqfilter.gaf > gapped-unitig-unrolled-hifi-resolved.gfa
awk '{if ($2 >= 100000) {sum += $2*$3; count += $2;}}END{print sum/count;}' < nodecovs-ont.csv
scripts/estimate_unique_local.py gapped-unitig-unrolled-hifi-resolved.gfa alns-ont-filter-trim.gaf 100000 30 0.8 > unique_nodes_ont_coverage.txt
# scripts/translate_uniques.py normal-hifi_connected_twice.gfa < unique_nodes_hifi.txt > translated_uniques.txt
# scripts/translate_nodes_by_seq.py normal-hifi_connected_twice.gfa unitig-unrolled-hifi-resolved.gfa < translated_uniques.txt > unique_nodes_ont_translated.txt
# cat unique_nodes_ont_coverage.txt unique_nodes_ont_translated.txt | sort | uniq > unique_nodes_ont.txt
scripts/fix_diploid_unique_nodes.py unique_nodes_ont_coverage.txt nodecovs-ont.csv gapped-unitig-unrolled-hifi-resolved.gfa > unique_nodes_diploidfix.txt
cp unique_nodes_diploidfix.txt unique_nodes_ont.txt

scripts/find_bridges.py unique_nodes_ont.txt < paths.txt > bridges.txt
grep -v '(' < bridges.txt | grep -vP '^$' | scripts/remove_wrong_connections_2.py | sort > bridging_seq_all.txt
scripts/pick_majority_bridge.py < bridging_seq_all.txt > bridging_seq_picked_all.txt
scripts/remove_crosslink_paths.py unique_nodes_ont.txt bridging_seq_picked_all.txt bridges.txt > bridges_fixcrosslink.txt 2> forbidden_crosslinks.txt
scripts/fix_diploid_paths.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridges_fixcrosslink.txt bridges.txt 3 > bridging_seq_diploidfix_all.txt
cp bridging_seq_diploidfix_all.txt bridging_seq_picked.txt
scripts/forbid_unbridged_tangles.py unique_nodes_ont.txt gapped-unitig-unrolled-hifi-resolved.gfa bridging_seq_all.txt paths.txt nodecovs-ont.csv 30 > forbidden_ends.txt
scripts/connect_uniques.py gapped-unitig-unrolled-hifi-resolved.gfa forbidden_ends.txt bridging_seq_picked.txt > connected.gfa

scripts/merge_unresolved_dbg_nodes.py < connected.gfa > normal-connected.gfa
scripts/add_fake_alignments.py unitig-unrolled-hifi-resolved.gfa normal-connected.gfa alns-ont-filter-trim.gaf nodecovs-ont.csv fake-ont-alns.gaf fake-ont-nodecovs.csv 10
/usr/bin/time -v scripts/resolve_triplets_kmerify.py normal-connected.gfa fake-ont-paths.txt fake-ont-nodecovs.csv 100000 3 5 3 2 < fake-ont-alns.gaf > ont-resolved-graph.gfa 2> stderr_ont_resolved_graph.txt

scripts/unroll_tip_loops.py ont-resolved-graph.gfa 3 < fake-ont-paths.txt > unrolled-ont-resolved.gfa
scripts/unitigify.py < unrolled-ont-resolved.gfa > unitig-unrolled-ont-resolved.gfa

# consensus
grep -P '^S' < unitig-unrolled-ont-resolved.gfa | awk '{print ">" $2; print $3;}' > contigs_rle.fa
winnowmap -x map-pb -t 32 contigs_rle.fa hifi.fa > alns.paf
scripts/pick_reads_stdin.py used_ont.txt < ont.fa > ont_gap_subset.fa
scripts/rle.py < ont_gap_subset.fa > ont_gap_subset_rle.fa
winnowmap -x map-ont -t 32 contigs_rle.fa ont_gap_subset_rle.fa >> alns.paf
scripts/get_layout_from_aln.py contigs_rle.fa alns.paf read_names.txt hifi.fa ont_gap_subset.fa > layout.txt

# just for debug info
scripts/check_layout_gaps.py < layout.txt > gaps.txt

# only for evaluating CHM13 haploid assemblies
# UntipRelative 30000 30000 0.1 0.1 < unitig-unrolled-ont-resolved.gfa > connected-tip.gfa
# scripts/unitigify.py < connected-tip.gfa > unitig-normal-connected-tip.gfa
# scripts/pop_bubbles_keep_longest.py 10 < unitig-normal-connected-tip.gfa > popped-unitig-normal-connected-tip.gfa
# scripts/unitigify.py < popped-unitig-normal-connected-tip.gfa > unitig-popped-unitig-normal-connected-tip.gfa
