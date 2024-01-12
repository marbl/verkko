#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
unique_nodes_file = sys.argv[2]
picked_paths_file = sys.argv[3]
counted_paths_file = sys.argv[4]
# paths to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:]

def startnode(p):
	path = p.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	return path[0]

def endnode(p):
	path = p.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	return revnode(path[-1])

def canontip(n1, n2):
	if n2 + n1 < n1 + n2: return (n2, n1)
	return (n1, n2)

def check_side(start, checked_sides, graph_edges):
	global forbidden_connections
	if start in checked_sides: return
	if start not in connections_per_unique: return
	tangle_sides = set()
	check_stack = [start]
	tangle_nodes = set()
	while len(check_stack) > 0:
		top = check_stack.pop()
		if top in tangle_nodes: continue
		tangle_nodes.add(top)
		if top[1:] in unique_nodes:
			tangle_sides.add(top)
		else:
			check_stack.append(revnode(top))
		if top in graph_edges:
			for edge in graph_edges[top]:
				check_stack.append(edge)
	for node in tangle_sides:
		assert node not in checked_sides
		checked_sides.add(node)
	if len(tangle_sides) % 2 != 0: return
	best_per_side = {}
	for node in tangle_sides:
		best_count = 0
		best_connection = None
		if node not in connections_per_unique: return
		for connection in connections_per_unique[node]:
			key = canontip(node, connection)
			if connection_counts[key] > best_count:
				best_count = connection_counts[key]
				best_connection = connection
			elif connection_counts[key] == best_count:
				best_connection = None
		if best_connection is None: return
		best_per_side[node] = best_connection
	for node in best_per_side:
		if best_per_side[best_per_side[node]] != node: return
	for node in tangle_sides:
		for connection in connections_per_unique[node]:
			if connection != best_per_side[node]:
				key = canontip(node, connection)
				forbidden_connections.add(key)

graph_edges = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			tonode = revnode(tonode)
			if fromnode not in graph_edges: graph_edges[fromnode] = set()
			if tonode not in graph_edges: graph_edges[tonode] = set()
			graph_edges[fromnode].add(tonode)
			graph_edges[tonode].add(fromnode)

unique_nodes = set()
with open(unique_nodes_file) as f:
	for l in f:
		unique_nodes.add(l.strip())

paths_per_connection = {}
connections_per_unique = {}
all_paths = set()
with open(picked_paths_file) as f:
	for l in f:
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		if len(path) == 0: continue
		all_paths.add(l.strip())
		if path[0] not in connections_per_unique: connections_per_unique[path[0]] = set()
		connections_per_unique[path[0]].add(revnode(path[-1]))
		if revnode(path[-1]) not in connections_per_unique: connections_per_unique[revnode(path[-1])] = set()
		connections_per_unique[revnode(path[-1])].add(path[0])
		key = canontip(path[0], path[-1])
		if key not in paths_per_connection: paths_per_connection[key] = set()
		paths_per_connection[key].add(l.strip())

connection_counts = {}
with open(counted_paths_file) as f:
	for l in f:
		if len(l.strip()) == 0: continue
		path = l.strip().replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		if len(path) == 0: continue
		key = canontip(path[0], revnode(path[-1]))
		if key not in connection_counts: connection_counts[key] = 0
		connection_counts[key] += 1

forbidden_connections = set()
checked_sides = set()
for node in unique_nodes:
	check_side('>' + node, checked_sides, graph_edges)
	check_side('<' + node, checked_sides, graph_edges)

for pathstr in all_paths:
	path = pathstr.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	if len(path) == 0: continue
	key = canontip(path[0], revnode(path[-1]))
	if key not in forbidden_connections:
		print(pathstr)
	else:
		sys.stderr.write(pathstr + "\n")
