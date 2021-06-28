#!/usr/bin/python

import sys

graph_file = sys.argv[1]
alignment_file = sys.argv[2]
node_coverage_file = sys.argv[3]
max_removable_length = int(sys.argv[4])
min_safe_coverage = float(sys.argv[5])
safe_path_coverage = int(sys.argv[6])
min_path_overlap_len = int(sys.argv[7])
# gfa graph to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(left, right):
	fwstr = left + right
	bwstr = revnode(right) + revnode(left)
	if bwstr < fwstr: return (revnode(right), revnode(left))
	return (left, right)

def path_overlap_len(path1, path2, offset, node_lens, edge_overlaps):
	matching_path = []
	for i in range(0, len(path1)+len(path2)):
		path1_pos = offset + i
		path2_pos = i
		if path1_pos < 0: continue
		if path1_pos >= len(path1): continue
		if path2_pos < 0: continue
		if path2_pos >= len(path2): continue
		if path1[path1_pos] != path2[path2_pos]: return 0
		matching_path.append(path1[path1_pos])
	if len(matching_path) == 0: return 0
	match_length = node_lens[matching_path[0][1:]]
	for i in range(1, len(matching_path)):
		overlap = edge_overlaps[canon(matching_path[i-1], matching_path[i])]
		match_length += node_lens[matching_path[i][1:]] - overlap
	return match_length

def paths_overlap(path1, path2, node_lens, edge_overlaps, min_path_overlap_len):
	max_path_overlap = 0
	rev_path2 = [revnode(n) for n in path2[::-1]]
	for i in range(-len(path2), len(path1)):
		overlap_here = path_overlap_len(path1, path2, i, node_lens, edge_overlaps)
		max_path_overlap = max(max_path_overlap, overlap_here)
		overlap_here = path_overlap_len(path1, rev_path2, i, node_lens, edge_overlaps)
		max_path_overlap = max(max_path_overlap, overlap_here)
	if max_path_overlap >= min_path_overlap_len: return True
	return False

node_lines = {}
node_lens = {}
edges = {}
edge_overlaps = {}
maybe_removable_node = set()

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if len(parts[2]) <= max_removable_length: maybe_removable_node.add(parts[1])
			node_lens[parts[1]] = len(parts[2])
			node_lines[parts[1]] = l.strip()
		if parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))
			edge_overlaps[canon(fromnode, tonode)] = int(parts[5][:-1])

paths_crossing = {}
with open(alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		path = parts[5].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		part_paths = []
		last_break = 0
		for i in range(0, len(path)):
			if path[i][1:] not in node_lens:
				if i > last_break: part_paths.append(path[last_break:i])
				last_break = i+1
			elif i > 0 and (path[i-1] not in edges or path[i] not in edges[path[i-1]]):
				if i > last_break: part_paths.append(path[last_break:i])
				last_break = i
		if last_break < len(path): part_paths.append(path[last_break:])
		for part in part_paths:
			path_nodes = set(n[1:] for n in path)
			for node in path_nodes:
				if node not in paths_crossing: paths_crossing[node] = []
				paths_crossing[node].append(part)

node_coverages = {}
with open(node_coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node": continue
		node = parts[0]
		coverage = float(parts[2])
		if coverage >= min_safe_coverage and node in maybe_removable_node: maybe_removable_node.remove(node)
		node_coverages[node] = coverage

removed_nodes = set()
for node in maybe_removable_node:
	if node in node_coverages and node_coverages[node] >= min_safe_coverage: continue
	has_safe_cluster = False
	if node in paths_crossing:
		parent = []
		cluster_count = []
		for i in range(0, len(paths_crossing[node])):
			parent.append(i)
			cluster_count.append(0)
		for i in range(0, len(paths_crossing[node])):
			for j in range(i+1, len(paths_crossing[node])):
				if paths_overlap(paths_crossing[node][i], paths_crossing[node][j], node_lens, edge_overlaps, min_path_overlap_len):
					parent[j] = i
		for i in range(0, len(paths_crossing[node])):
			while parent[i] != parent[parent[i]]: parent[i] = parent[parent[i]]
		for i in range(0, len(paths_crossing[node])):
			cluster_count[parent[i]] += 1
		for i in range(0, len(paths_crossing[node])):
			if cluster_count[i] >= safe_path_coverage: has_safe_cluster = True
	if has_safe_cluster: continue
	removed_nodes.add(node)

sys.stderr.write("removed " + str(len(removed_nodes)) + " nodes" + "\n")

for node in node_lines:
	if node in removed_nodes: continue
	print(node_lines[node])

for edge in edges:
	if edge[1:] in removed_nodes: continue
	for target in edges[edge]:
		if target[1:] in removed_nodes: continue
		overlap = edge_overlaps[canon(edge, target)]
		print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + target[1:] + "\t" + ("+" if target[0] == ">" else "-") + "\t" + str(overlap) + "M")

