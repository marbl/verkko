#!/usr/bin/env python

import sys
import graph_functions as gf

forbidden_ends_file = sys.argv[1]
bridges_file = sys.argv[2]
old_coverage_csv_in = sys.argv[3]
new_coverage_csv_out = sys.argv[4]
extra_coverage = int(sys.argv[5])
# fake gaf to stdout

forbidden_ends = set()
with open(forbidden_ends_file) as f:
	for l in f:
		forbidden_ends.add(l.strip())

nodes_with_extra_coverage = set()
with open(bridges_file) as f:
	for l in f:
		path = l.strip().replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		if path[0] in forbidden_ends or gf.revnode(path[-1]) in forbidden_ends: continue
		for node in path:
			nodes_with_extra_coverage.add(node[1:])
		for i in range(0, extra_coverage):
			print("FAKE\t1\t0\t1\t+\t" + "".join(path) + "\t1\t0\t1\t1\t1\t60\tNM:i:0\tAS:f:100000\tdv:f:0\tid:f:1")

with open(old_coverage_csv_in) as f:
	with open(new_coverage_csv_out, "w") as out:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] in nodes_with_extra_coverage:
				parts[2] = str(float(parts[2]) + extra_coverage)
			out.write("\t".join(parts) + "\n")


