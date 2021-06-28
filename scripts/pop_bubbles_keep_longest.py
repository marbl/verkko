#!/usr/bin/python

import sys

max_bubble_size = int(sys.argv[1])
# gfa from stdin
# gfa to stdout

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def getone(s):
	assert len(s) >= 1
	for n in s:
		return n

def canontip(left, right):
	fwstr = left + right
	bwstr = right + left
	if bwstr < fwstr: return (right, left)
	return (left, right)

def remove_graph_node(node, nodeseqs, edges):
	assert node in nodeseqs
	del nodeseqs[node]
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

def pop_bubble(start, end, nodeseqs, edges, edge_overlaps):
	bubble_nodes = set()
	stack = [(start, start, 0)]
	longest = {}
	while len(stack) > 0:
		(n, predecessor, length) = stack.pop()
		if n in longest and longest[n][0] >= length: continue
		longest[n] = (length, predecessor)
		if n != start: bubble_nodes.add(n)
		if n == end: continue
		assert n in edges
		assert len(edges[n]) >= 1
		for edge in edges[n]: stack.append((edge, n, length + len(nodeseqs[edge[1:]]) - edge_overlaps[canontip(n, revnode(edge))]))
	kept = [end]
	while kept[-1] != start:
		kept.append(longest[kept[-1]][1])
	kept = kept[::-1]
	assert kept[0] == start
	assert kept[-1] == end
	kept_set = set(kept)
	for node in bubble_nodes:
		if node in kept_set: continue
		remove_graph_node(node[1:], nodeseqs, edges)
	for i in range(0, len(kept)):
		if i != 0:
			if len(edges[revnode(kept[i])]) >= 2:
				sys.stderr.write(kept[i] + "\n")
			del edges[revnode(kept[i])]
		if i != len(kept)-1:
			if len(edges[kept[i]]) >= 2:
				sys.stderr.write(kept[i] + "\n")
			del edges[kept[i]]
	for i in range(1, len(kept)):
		old = kept[i-1]
		current = kept[i]
		edges[old] = set()
		edges[old].add(current)
		edges[revnode(current)] = set()
		edges[revnode(current)].add(revnode(old))

nodeseqs = {}
edges = {}
edge_overlaps = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == 'S':
		nodeseqs[parts[1]] = parts[2]
	elif parts[0] == 'L':
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		edges[fromnode].add(tonode)
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[revnode(tonode)].add(revnode(fromnode))
		edge_overlaps[canontip(fromnode, revnode(tonode))] = int(parts[5][:-1])

check_nodes = set(nodeseqs)

for node in check_nodes:
	if node not in nodeseqs: continue
	bubble = find_bubble(">" + node, edges)
	if bubble is not None:
		pop_bubble(bubble[0], bubble[1], nodeseqs, edges, edge_overlaps)
	bubble = find_bubble("<" + node, edges)
	if bubble is not None:
		pop_bubble(bubble[0], bubble[1], nodeseqs, edges, edge_overlaps)

for n in nodeseqs:
	print("S\t" + n + "\t" + nodeseqs[n])
for edge in edges:
	for target in edges[edge]:
		key = canontip(edge, revnode(target))
		overlap = edge_overlaps[key]
		assert key in edge_overlaps
		print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + target[1:] + "\t" + ("+" if target[0] == ">" else "-") + "\t" + str(overlap) + "M")
