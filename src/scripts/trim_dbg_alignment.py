#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
edge_trim = 0
if len(sys.argv) >= 3: edge_trim = int(sys.argv[2])
# gfa from stdin
# gfa to stdout

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

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
			edge_overlaps[(revnode(tonode), revnode(fromnode))] = int(parts[5][:-1])

for l in sys.stdin:
	parts = l.strip().split('\t')
	path = parts[5].replace("<", "\t<").replace(">", "\t>").strip().split("\t")
	assert len(path) >= 1
	parts[2] = str(int(parts[2]) + edge_trim)
	parts[3] = str(int(parts[3]) - edge_trim)
	if int(parts[3]) <= int(parts[2]): continue
	parts[7] = str(int(parts[7]) + edge_trim)
	parts[8] = str(int(parts[8]) - edge_trim)
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
