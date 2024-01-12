#!/usr/bin/env python

import sys

original_selected = sys.argv[1]
gap_paths = sys.argv[2]
alignment_file = sys.argv[3]
orig_alignment_file = sys.argv[4]
min_q = int(sys.argv[5])

selected = set()
with open(original_selected) as f:
	for l in f:
		selected.add(l.strip())

# record reads which are used for gap patching
gap_reads = set()
with open(gap_paths) as f:
	for l in f:
		parts = l.strip().split('\t')
		gap_reads.add(parts[0])

# record which nodes were cut
cut_nodes = set()
with open(alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		for p in parts[5].replace(">", " ").replace("<", " ").strip().split(' '):
			if "cut" in p:
				cut_nodes.add(p.split("_")[0])
				gap_reads.add(parts[0])

# now go through and add any mappings that aren't using cut nodes or reads in a gap
with open(orig_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		nodes = parts[5].replace(">", " ").replace("<", " ").strip().split(' ')
		if int(parts[11]) >= min_q and parts[0] in selected and parts[0] not in gap_reads and len(cut_nodes.intersection(nodes)) == 0:
			print(l.strip())
