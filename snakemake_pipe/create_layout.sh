#!/usr/bin/sh

grep -P '^S' < unitig-popped-unitig-normal-connected-tip.gfa | awk '{print ">" $2; print $3;}' > contigs_rle.fa
/usr/bin/time -v /data/Phillippy/tools/Winnowmap/bin/winnowmap -x map-pb -t 32 contigs_rle.fa ../chm13_consensus_test/m*.fasta > alns.paf 2> stderr_winnowmap_hifi.txt
zcat /data/Phillippy/seq/chm13/nanopore/split/*.fq.gz | awk '{if (NR % 4 == 1 || NR % 4 == 2) print;}' | tr '@' '>' | /data/rautiainenma/hybrid-assembly/scripts/pick_reads_stdin.py used_ont.txt > ont_gap_subset.fa
/data/rautiainenma/hybrid-assembly/scripts/rle.py < ont_gap_subset.fa > ont_gap_subset_rle.fa
/usr/bin/time -v /data/Phillippy/tools/Winnowmap/bin/winnowmap -x map-ont -t 32 contigs_rle.fa ont_gap_subset_rle.fa >> alns.paf 2> stderr_winnowmap_ont.txt
/usr/bin/time -v /data/rautiainenma/hybrid-assembly/scripts/get_layout_from_aln.py contigs_rle.fa alns.paf read_names.txt ../chm13_consensus_test/m*.fasta ont_gap_subset.fa > layout.txt 2> stderr_layout.txt

/data/rautiainenma/hybrid-assembly/scripts/check_layout_gaps.py < layout.txt > gaps.txt

