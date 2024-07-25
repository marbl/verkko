#!/usr/bin/env python

import sys
import graph_functions as gf

unique_nodes_file = sys.argv[1]
graph_file = sys.argv[2]
picked_paths_file = sys.argv[3]
counted_paths_file = sys.argv[4]
other_path_create_coverage_threshold = int(sys.argv[5])
# new picked paths to stdout

unique_nodes = set()
with open(unique_nodes_file) as f:
	for l in f:
		unique_nodes.add(l.strip())

nodes = set()
edges = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodes.add(parts[1])
		elif parts[0] == 'L':
			fromnode = ('>' if parts[2] == '+' else '<') + parts[1]
			tonode = ('>' if parts[4] == '+' else '<') + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
			edges[gf.revnode(tonode)].add(gf.revnode(fromnode))

check_nodes = set()
discarded_paths_per_checkable = {}
for node in nodes:
	if node in unique_nodes: continue
	if '>' + node not in edges: continue
	if '<' + node not in edges: continue
	if len(edges['>' + node]) != 2: continue
	if len(edges['<' + node]) != 2: continue
	diploid_core = True
	neighbors = set()
	for edge in edges['>' + node]:
		if edge[1:] not in unique_nodes:
			diploid_core = False
			break
		neighbors.add(edge[1:])
	for edge in edges['<' + node]:
		if edge[1:] not in unique_nodes:
			diploid_core = False
			break
		neighbors.add(edge[1:])
	if not diploid_core: continue
	if len(neighbors) != 4: continue
	check_nodes.add(node)
	discarded_paths_per_checkable[node] = set()

paths_per_checkable = {}
checkable_paths = set()
with open(picked_paths_file) as f:
	for l in f:
		parts = l.replace('>', '\t').replace('<', '\t').strip().split('\t')
		is_checkable = False
		for node in parts:
			if node not in check_nodes: continue
			assert len(parts) == 3
			assert parts[0] != node
			assert parts[1] == node
			assert parts[2] != node
			checkable_paths.add(l.strip())
			if node not in paths_per_checkable: paths_per_checkable[node] = set()
			paths_per_checkable[node].add(l.strip())
			is_checkable = True
		if not is_checkable: print(l.strip())

path_counts = {}
with open(counted_paths_file) as f:
	for l in f:
		is_checkable = False
		for node in l.replace('>', '\t').replace('<', '\t').strip().split('\t'):
			if node not in check_nodes: continue
			is_checkable = True
		if not is_checkable: continue
		if l.strip() not in path_counts: path_counts[l.strip()] = 0
		path_counts[l.strip()] += 1
		if l.strip() not in checkable_paths:
			for node in l.replace('>', '\t').replace('<', '\t').strip().split('\t'):
				if node not in check_nodes: continue
				discarded_paths_per_checkable[node].add(l.strip())

for node in check_nodes:
	if node not in paths_per_checkable: continue
	if len(paths_per_checkable[node]) == 1:
		if path_counts[gf.getone(paths_per_checkable[node])] >= other_path_create_coverage_threshold and len(discarded_paths_per_checkable[node]) == 0:
			used_nodes = set()
			for part in gf.getone(paths_per_checkable[node]).replace('>', '\t').replace('<', '\t').strip().split('\t'):
				used_nodes.add(part)
			assert len(used_nodes) == 3
			assert node in used_nodes
			prev_node = None
			next_node = None
			for edge in edges['<' + node]:
				if edge[1:] in used_nodes: continue
				assert prev_node is None
				prev_node = gf.revnode(edge)
			for edge in edges['>' + node]:
				if edge[1:] in used_nodes: continue
				assert next_node is None
				next_node = edge
			assert prev_node is not None
			assert next_node is not None
			print(prev_node + '>' + node + next_node)
		for path in paths_per_checkable[node]: print(path)
		continue
	elif len(paths_per_checkable[node]) == 2:
		for path in paths_per_checkable[node]: print(path)
		continue
	elif len(paths_per_checkable[node]) >= 3:
		paths_and_covs = []
		for path in paths_per_checkable[node]:
			assert path in path_counts
			paths_and_covs.append((path, path_counts[path]))
		paths_and_covs.sort(key=lambda x: -x[1])
		while len(paths_and_covs) > 2 and paths_and_covs[-1][1] < paths_and_covs[1][1]: paths_and_covs.pop()
		valid_options = []
		assert len(paths_and_covs) >= 2
		for i in range(0, len(paths_and_covs)):
			if paths_and_covs[i][1] < paths_and_covs[0][1]: break
			for j in range(i+1, len(paths_and_covs)):
				path1 = paths_and_covs[i][0].replace('>', '\t').replace('<', '\t').strip().split('\t')
				path2 = paths_and_covs[j][0].replace('>', '\t').replace('<', '\t').strip().split('\t')
				assert len(path1) == 3
				assert len(path2) == 3
				used_nodes = set()
				for used in path1: used_nodes.add(used)
				for used in path2: used_nodes.add(used)
				assert len(used_nodes) <= 5
				if len(used_nodes) == 5: valid_options.append((paths_and_covs[i][0], paths_and_covs[j][0]))
		if len(valid_options) == 1:
			print(valid_options[0][0])
			print(valid_options[0][1])
		else:
			for path in paths_per_checkable[node]: print(path)

