#!/usr/bin/env python

import sys
import graph_functions as gf

graph_file = sys.argv[1]
edge_trim = 0
ignore_file = ""
if len(sys.argv) >= 3: edge_trim = int(sys.argv[2])
if len(sys.argv) >= 4: ignore_file = sys.argv[3] 
# gaf from stdin
# gaf to stdout

ignore_reads = set()
if ignore_file != "":
	with open(ignore_file) as f:
		for l in f:
			ignore_reads.add(l.strip().split('\t')[0])

edge_overlaps = {}
node_lens = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			node_lens[parts[1]] = len(parts[2])
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			edge_overlaps[(fromnode, tonode)] = int(parts[5][:-1])
			edge_overlaps[(gf.revnode(tonode), gf.revnode(fromnode))] = int(parts[5][:-1])

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] in ignore_reads: continue
	path = parts[5].replace("<", "\t<").replace(">", "\t>").strip().split("\t")
	assert len(path) >= 1
	while len(path) >= 2 and int(parts[7]) >= node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]:
		parts[6] = str(int(parts[6]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		parts[7] = str(int(parts[7]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		parts[8] = str(int(parts[8]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		path = path[1:]
	while len(path) >= 2 and int(parts[6]) - int(parts[8]) >= node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])]:
		parts[6] = str(int(parts[6]) - (node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])]))
		path = path[:-1]
	back_trim = 0
	front_trim = 0
	if len(path) >= 2 and int(parts[7]) + edge_trim >= node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]:
		back_trim = node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])] - int(parts[7])
	if len(path) >= 2 and int(parts[6]) - int(parts[8]) + edge_trim >= node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])]:
		front_trim = node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])] - (int(parts[6]) - int(parts[8]))
	parts[2] = str(int(parts[2]) + back_trim)
	parts[3] = str(int(parts[3]) - front_trim)
	if int(parts[3]) <= int(parts[2]): continue
	parts[7] = str(int(parts[7]) + back_trim)
	parts[8] = str(int(parts[8]) - front_trim)
	while len(path) >= 2 and int(parts[7]) >= node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]:
		parts[6] = str(int(parts[6]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		parts[7] = str(int(parts[7]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		parts[8] = str(int(parts[8]) - (node_lens[path[0][1:]] - edge_overlaps[(path[0], path[1])]))
		path = path[1:]
	while len(path) >= 2 and int(parts[6]) - int(parts[8]) >= node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])]:
		parts[6] = str(int(parts[6]) - (node_lens[path[-1][1:]] - edge_overlaps[(path[-2], path[-1])]))
		path = path[:-1]
	assert len(path) >= 1
	parts[5] = "".join(path)
	print("\t".join(parts))
