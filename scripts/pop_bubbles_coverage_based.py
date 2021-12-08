#!/usr/bin/env python

import sys

node_coverage_file = sys.argv[1]
# gfa from stdin
# gfa to stdout

max_bubble_size = 10

def getone(s):
	for n in s:
		return n
	assert False

def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

def remove_graph_node(node, edges):
	if ">" + node in edges:
		for edge in edges[">" + node]:
			assert revnode(edge) in edges
			assert "<" + node in edges[revnode(edge)]
			edges[revnode(edge)].remove("<" + node)
		del edges[">" + node]
	if "<" + node in edges:
		for edge in edges["<" + node]:
			assert revnode(edge) in edges
			assert ">" + node in edges[revnode(edge)]
			edges[revnode(edge)].remove(">" + node)
		del edges["<" + node]

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013
# fig. 5
def find_bubble(s, edges):
	if s not in edges: return None
	if len(edges[s]) < 2: return None
	S = [s]
	visited = set()
	seen = set()
	seen.add(s)
	while len(S) > 0:
		v = S.pop()
		assert v in seen
		seen.remove(v)
		assert v not in visited
		visited.add(v)
		if len(visited) > max_bubble_size: return None
		if v not in edges: return None
		if len(edges[v]) == 0: return None
		for u in edges[v]:
			if u[1:] == v[1:]: return None
			if revnode(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert revnode(u) in edges
			assert len(edges[revnode(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in edges[revnode(u)]:
				parent = revnode(parent_edge)
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == getone(seen):
			t = S.pop()
			if t in edges:
				for edge in edges[t]:
					if edge == s: return None
			return (s, t)
	return None

def pop_bubble(start, end, removed_nodes, removed_edges, edges, coverage):
	bubble_nodes = set()
	bubble_edges = set()
	predecessor = {}
	stack = []
	stack.append((start, start))
	while len(stack) > 0:
		(top, before) = stack.pop()
		if top not in predecessor: predecessor[top] = before
		if before in coverage and predecessor[top] not in coverage: predecessor[top] = before
		if before in coverage and predecessor[top] in coverage and coverage[before] > coverage[predecessor[top]]: predecessor[top] = before
		bubble_nodes.add(top[1:])
		bubble_edges.add((before, top))
		if top == end: continue
		for edge in edges[top]:
			stack.append((edge, top))
	assert end in predecessor
	path = [end]
	while path[-1] != start:
		path.append(predecessor[path[-1]])
	path = path[::-1]
	kept_nodes = set()
	kept_edges = set()
	for node in path:
		kept_nodes.add(node[1:])
	for i in range(1, len(path)):
		kept_edges.add((path[i-1], path[i]))
	assert len(kept_nodes) == len(path)
	assert len(kept_edges) == len(path)-1
	for node in bubble_nodes:
		if node in kept_nodes: continue
		remove_graph_node(node, edges)
		removed_nodes.add(node)
	for edge in bubble_edges:
		if edge in kept_edges: continue
		if (revnode(edge[1]), revnode(edge[0])) in kept_edges: continue
		removed_edges.add(edge)
		if edge[0] in edges:
			if edge[1] in edges[edge[0]]:
				edges[edge[0]].remove(edge[1])
		if revnode(edge[1]) in edges:
			if revnode(edge[0]) in edges[revnode(edge[1])]:
				edges[revnode(edge[1])].remove(revnode(edge[0]))

coverage = {}
with open(node_coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node" and parts[1] == "coverage": continue
		coverage[parts[0]] = float(parts[1])

nodelens = {}
edges = {}
nodelines = []
edgelines = []
removed_nodes = set()
removed_edges = set()

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == 'S':
		nodelines.append((parts[1], l.strip()))
		nodelens[parts[1]] = len(parts[2])
	elif parts[0] == 'L':
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		edgelines.append((fromnode, tonode, l.strip()))
		if fromnode not in edges: edges[fromnode] = set()
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[fromnode].add(tonode)
		edges[revnode(tonode)].add(revnode(fromnode))

long_coverage_cov_sum = 0.0
long_coverage_len_sum = 0.0
for node in coverage:
	assert node in nodelens
	if nodelens[node] < 100000: continue
	long_coverage_len_sum += nodelens[node]
	long_coverage_cov_sum += nodelens[node] * coverage[node]

avg_coverage = long_coverage_cov_sum / long_coverage_len_sum
sys.stderr.write("average coverage " + str(avg_coverage) + "\n")

for node in coverage:
	if node in removed_nodes: continue
	if coverage[node] < avg_coverage * 0.5: continue
	if coverage[node] > avg_coverage * 1.5: continue
	bubble = find_bubble(">" + node, edges)
	if bubble:
		assert bubble[0] == ">" + node
		assert bubble[1][1:] != node
		if bubble[1][1:] in coverage and coverage[bubble[1][1:]] >= avg_coverage * 0.5 and coverage[bubble[1][1:]] <= avg_coverage * 1.5:
			pop_bubble(bubble[0], bubble[1], removed_nodes, removed_edges, edges, coverage)
	bubble = find_bubble("<" + node, edges)
	if bubble:
		assert bubble[0] == "<" + node
		assert bubble[1][1:] != node
		if bubble[1][1:] in coverage and coverage[bubble[1][1:]] >= avg_coverage * 0.5 and coverage[bubble[1][1:]] <= avg_coverage * 1.5:
			pop_bubble(bubble[0], bubble[1], removed_nodes, removed_edges, edges, coverage)

for node in removed_nodes:
	sys.stderr.write(node + "\n")

for edge in removed_edges:
	if edge[0][1:] in removed_nodes: continue
	if edge[1][1:] in removed_nodes: continue
	sys.stderr.write(edge[0] + "\t" + edge[1] + "\n")

for node in nodelines:
	if node[0] in removed_nodes: continue
	print(node[1])

for edge in edgelines:
	if edge[0][1:] in removed_nodes: continue
	if edge[1][1:] in removed_nodes: continue
	if (edge[0], edge[1]) in removed_edges: continue
	if (revnode(edge[1]), revnode(edge[0])) in removed_edges: continue
	print(edge[2])
