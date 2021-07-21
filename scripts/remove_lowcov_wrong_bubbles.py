#!/usr/bin/python

import sys

graph_file = sys.argv[1]
coverage_file = sys.argv[2]
# alignment_file = sys.argv[3]
max_removable_coverage = float(sys.argv[3])
max_removable_len = int(sys.argv[4])
min_safe_coverage = int(sys.argv[5])
# gfa to stdout

possibly_removable_nodes = set()
safe_nodes = set()
with open(coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node": continue
		if int(parts[1]) <= max_removable_len and float(parts[2]) <= max_removable_coverage: possibly_removable_nodes.add(parts[0])
		if float(parts[2]) >= min_safe_coverage: safe_nodes.add(parts[0])

node_lines = []
edge_lines = []
edges = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			node_lines.append((parts[1], l.strip()))
		if parts[0] == 'L':
			edge_lines.append((parts[1], parts[3], l.strip()))
			if (parts[1], parts[2]) not in edges: edges[(parts[1], parts[2])] = set()
			edges[(parts[1], parts[2])].add((parts[3], parts[4]))
			if (parts[3], ("+" if parts[4] == "-" else "-")) not in edges: edges[(parts[3], ("+" if parts[4] == "-" else "-"))] = set()
			edges[(parts[3], ("+" if parts[4] == "-" else "-"))].add((parts[1], ("+" if parts[2] == "-" else "-")))

removed = set()
for node in possibly_removable_nodes:
	if (node, "+") not in edges: continue
	if (node, "-") not in edges: continue
	if len(edges[(node, "+")]) >= 2: continue
	if len(edges[(node, "-")]) >= 2: continue
	removable = True
	if (node, "+") in edges:
		for edge in edges[(node, "+")]:
			if edge[0] not in safe_nodes: removable = False
			revkey = (edge[0], "+" if edge[1] == "-" else "-")
			assert revkey in edges
			has_other_solid = False
			for edge2 in edges[revkey]:
				if edge2[0] in safe_nodes: has_other_solid = True
			if not has_other_solid: removable = False
	if (node, "-") in edges:
		for edge in edges[(node, "-")]:
			if edge[0] not in safe_nodes: removable = False
			revkey = (edge[0], "+" if edge[1] == "-" else "-")
			assert revkey in edges
			has_other_solid = False
			for edge2 in edges[revkey]:
				if edge2[0] in safe_nodes: has_other_solid = True
			if not has_other_solid: removable = False
	if not removable: continue
	removed.add(node)
	sys.stderr.write(node + "\n")

for n in node_lines:
	if n[0] in removed: continue
	print(n[1])

for e in edge_lines:
	if e[0] in removed: continue
	if e[1] in removed: continue
	print(e[2])
