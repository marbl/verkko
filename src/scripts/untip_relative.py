#!/usr/bin/env python

import sys
import math

max_removable_len = int(sys.argv[1])
min_safe_len = int(sys.argv[2])
fraction = float(sys.argv[3])
# input gfa from stdin
# output gfa to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

# https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
# iterative so we don't hit max recursion depth
def strong_connect_iterative(start, i, index, lowlink, on_stack, result, edges):
	stack = []
	stack.append((0, start, 0))
	S = []
	while len(stack) > 0:
		(state, node, neighbor_i) = stack[-1]
		stack.pop()
		if state == 0:
			assert node not in on_stack
			assert node not in index
			assert node not in lowlink
			index[node] = i
			lowlink[node] = i
			i += 1
			S.append(node)
			on_stack.add(node)
			stack.append((1, node, 0))
		elif state == 1:
			if node in edges and neighbor_i < len(edges[node]):
				neighbor = edges[node][neighbor_i][0]
				if neighbor not in index:
					stack.append((2, node, neighbor_i))
					stack.append((0, neighbor, 0))
					continue
				elif neighbor in on_stack:
					assert neighbor in index
					assert node in lowlink
					lowlink[node] = min(lowlink[node], index[neighbor])
				neighbor_i += 1
			if node in edges and neighbor_i < len(edges[node]):
				stack.append((1, node, neighbor_i))
			else:
				stack.append((3, node, 0))
		elif state == 2:
			neighbor = edges[node][neighbor_i][0]
			assert neighbor in lowlink
			lowlink[node] = min(lowlink[node], lowlink[neighbor])
			neighbor_i += 1
			stack.append((1, node, neighbor_i))
		elif state == 3:
			assert node in lowlink
			assert node in index
			if lowlink[node] == index[node]:
				result.append([])
				stacknode = ""
				while stacknode != node:
					assert len(S) > 0
					stacknode = S[-1]
					S.pop()
					assert stacknode in on_stack
					on_stack.remove(stacknode)
					result[-1].append(stacknode)
		else:
			assert False
	assert len(S) == 0
	return i

# https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
def topological_sort(nodelens, edges):
	index = {}
	lowlink = {}
	on_stack = set()
	result = []
	i = 0
	for node in nodelens:
		if ">" + node not in index:
			old_i = i
			i = strong_connect_iterative(">" + node, i, index, lowlink, on_stack, result, edges)
			assert i > old_i
		if "<" + node not in index:
			old_i = i
			i = strong_connect_iterative("<" + node, i, index, lowlink, on_stack, result, edges)
			assert i > old_i
	result = result[::-1]
	belongs_to_component = {}
	for i in range(0, len(result)):
		assert len(result[i]) >= 1
		for node in result[i]:
			assert node not in belongs_to_component
			belongs_to_component[node] = i
	for node in nodelens:
		assert ">" + node in belongs_to_component
		assert "<" + node in belongs_to_component
	for node in belongs_to_component:
		if node in edges:
			for edge in edges[node]:
				assert belongs_to_component[edge[0]] >= belongs_to_component[node]
	return (result, belongs_to_component)

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
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[revnode(tonode)].add((revnode(fromnode), overlap))

# need a list because strong_connect_iterative starts and stops iterating in a non-standard order so sets won't work
for node in edges:
	edges[node] = list(edges[node])

(order, belongs_to_component) = topological_sort(nodelens, edges)
depths = get_node_depths(order, belongs_to_component, nodelens, edges)
kept = get_keepers(order, belongs_to_component, depths, nodelens, edges)

for line in nodelines:
	if line[0] in kept: print(line[1])

for line in edgelines:
	if line[0] in kept and line[1] in kept: print(line[2])
