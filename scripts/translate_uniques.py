#!/usr/bin/env python

import sys

new_graph_file = sys.argv[1]
# old uniques from stdin
# new uniques to stdout

old_uniques = set()
for l in sys.stdin:
	old_uniques.add(l.strip())

with open(new_graph_file) as f:
	for l in f:
		if l[0] == 'S':
			parts = l.strip().split('\t')
			node_nodes = parts[1].split('_')[3:]
			for node in node_nodes:
				if len(node) < 2: continue
				if node[-2:] != "nf" and node[-2:] != "nb": continue
				actual_node = node[:-2]
				if actual_node not in old_uniques: continue
				print(parts[1].strip())
				break

