#!/usr/bin/python

import sys

graph_file = sys.argv[1]
node_coverage_file = sys.argv[2]
gaf_file = sys.argv[3]
min_edge_coverage = int(sys.argv[4])
max_tip_length = int(sys.argv[5])
max_tip_coverage = int(sys.argv[6])
min_alt_coverage = int(sys.argv[7])
# graph to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(left, right):
	fwstr = left + right
	bwstr = revnode(right) + revnode(left)
	if bwstr < fwstr: return (revnode(right), revnode(left))
	return (left, right)

node_lines = {}
edges = {}
edge_overlaps = {}
maybe_removable_node = set()

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if len(parts[2]) <= max_tip_length: maybe_removable_node.add(parts[1])
			node_lines[parts[1]] = l.strip()
		if parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))
			edge_overlaps[canon(fromnode, tonode)] = parts[5]

edge_coverage = {}
with open(gaf_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		path = parts[5].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		for i in range(1, len(path)):
			key = canon(path[i-1], path[i])
			if key not in edge_coverage: edge_coverage[key] = 0
			edge_coverage[key] += 1

node_coverages = {}
with open(node_coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node": continue
		node = parts[0]
		coverage = float(parts[2])
		if coverage > max_tip_coverage and node in maybe_removable_node: maybe_removable_node.remove(node)
		node_coverages[node] = coverage

remove_edges = set()
for edge in edges:
	for target in edges[edge]:
		key = canon(edge, target)
		if key in edge_coverage and edge_coverage[key] >= min_edge_coverage: continue
		if node_coverages[key[0][1:]] < min_alt_coverage and node_coverages[key[1][1:]] < min_alt_coverage: continue
		remove_edges.add(key)

for key in remove_edges:
	edges[key[0]].remove(key[1])
	edges[revnode(key[1])].remove(revnode(key[0]))

removable_tips = set()

for node in maybe_removable_node:
	if ">" + node in edges and len(edges[">" + node]) > 0:
		if "<" + node in edges and len(edges["<" + node]) > 0:
			continue
	removable_tips.add(node)

removed_tips = set()

for node in removable_tips:
	assert ">" + node not in edges or "<" + node not in edges or len(edges[">" + node]) == 0 or len(edges["<" + node]) == 0
	has_valid_remover = False
	if ">" + node in edges and len(edges[">" + node]) > 0:
		for edge in edges[">" + node]:
			if node_coverages[edge[1:]] < min_alt_coverage: continue
			assert revnode(edge) in edges
			for edge2 in edges[revnode(edge)]:
				if edge2[1:] == node: continue
				if node_coverages[edge2[1:]] >= min_alt_coverage: has_valid_remover = True
	if "<" + node in edges and len(edges["<" + node]) > 0:
		for edge in edges["<" + node]:
			if node_coverages[edge[1:]] < min_alt_coverage: continue
			assert revnode(edge) in edges
			for edge2 in edges[revnode(edge)]:
				if edge2[1:] == node: continue
				if node_coverages[edge2[1:]] >= min_alt_coverage: has_valid_remover = True
	if has_valid_remover: removed_tips.add(node)

sys.stderr.write("removed " + str(len(removed_tips)) + " tips" + "\n")

for node in node_lines:
	if node in removed_tips: continue
	print(node_lines[node])

for edge in edges:
	if edge[1:] in removed_tips: continue
	for target in edges[edge]:
		if target[1:] in removed_tips: continue
		overlap = edge_overlaps[canon(edge, target)]
		print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + target[1:] + "\t" + ("+" if target[0] == ">" else "-") + "\t" + overlap)

