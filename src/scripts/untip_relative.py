#!/usr/bin/env python

import sys
import math
import graph_functions as gf

max_removable_len = int(sys.argv[1])
min_safe_len = int(sys.argv[2])
fraction = float(sys.argv[3])
# input gfa from stdin
# output gfa to stdout

def get_node_depths(order, belongs_to_component, nodelens, edges):
	depths = {}
	for i in range(len(order)-1, -1, -1):
		if len(order[i]) >= 2:
			for node in order[i]:
				depths[node] = math.inf
		else:
			assert len(order[i]) == 1
			node = order[i][0]
			depths[node] = 0
			if node in edges:
				for edge in edges[node]:
					assert belongs_to_component[edge[0]] >= belongs_to_component[node]
					if edge[0] == node:
						assert belongs_to_component[edge[0]] == belongs_to_component[node]
						depths[node] = math.inf
					else:
						assert belongs_to_component[edge[0]] > belongs_to_component[node]
						if depths[edge[0]] == math.inf:
							depths[node] = math.inf
						else:
							depths[node] = max(depths[node], depths[edge[0]] + nodelens[edge[0][1:]] - edge[1])
	return depths

def remove_rec(kept, node, edges):
	if node not in kept: return
	kept.remove(node)
	if node in edges:
		for edge in edges[node]:
			remove_rec(kept, edge[0], edges)

def get_keepers(order, belongs_to_component, depths, nodelens, edges):
	kept = set()
	for node in nodelens:
		kept.add(">" + node)
		kept.add("<" + node)
	for node in edges:
		max_out_length = 0
		for edge in edges[node]:
			depth_here = depths[edge[0]] + nodelens[edge[0][1:]] - edge[1]
			max_out_length = max(max_out_length, depth_here)
		if max_out_length < min_safe_len: continue
		for edge in edges[node]:
			depth_here = depths[edge[0]] + nodelens[edge[0][1:]] - edge[1]
			if depth_here >= max_removable_len: continue
			if depth_here > max_out_length * fraction: continue
			remove_rec(kept, edge[0], edges)
	really_kept = set()
	for node in nodelens:
		if ">" + node in kept and "<" + node in kept:
			really_kept.add(node)
	return really_kept


nodelines = []
edgelines = []
nodelens = {}
edges = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		nodelines.append((parts[1], l.strip()))
		nodelens[parts[1]] = len(parts[2])
	elif parts[0] == "L":
		edgelines.append((parts[1], parts[3], l.strip()))
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		assert parts[5][-1] == "M"
		overlap = int(parts[5][:-1])
		if fromnode not in edges: edges[fromnode] = set()
		edges[fromnode].add((tonode, overlap))
		if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
		edges[gf.revnode(tonode)].add((gf.revnode(fromnode), overlap))

# need a list because strong_connect_iterative starts and stops iterating in a non-standard order so sets won't work
for node in edges:
	edges[node] = list(edges[node])

(order, belongs_to_component) = gf.topological_sort(nodelens, edges)
depths = get_node_depths(order, belongs_to_component, nodelens, edges)
kept = get_keepers(order, belongs_to_component, depths, nodelens, edges)

for line in nodelines:
	if line[0] in kept: print(line[1])

for line in edgelines:
	if line[0] in kept and line[1] in kept: print(line[2])
