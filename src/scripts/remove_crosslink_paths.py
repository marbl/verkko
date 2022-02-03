#!/usr/bin/env python

import sys

unique_nodes_file = sys.argv[1]
picked_paths_file = sys.argv[2]
counted_paths_file = sys.argv[3]
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

def check_side(start):
	global forbidden_connections
	if start not in connections_per_unique: return
	if len(connections_per_unique[start]) != 2: return
	other_side = set()
	used_connections = set()
	for node in connections_per_unique[start]:
		other_side.add(node)
		used_connections.add(canontip(start, node))
	if len(other_side) != 2: return
	has_two = False
	this_side = set()
	has_many = False
	while True:
		added_any = False
		for side in other_side:
			for path in connections_per_unique[side]:
				if len(connections_per_unique[side]) == 2: has_two = True
				if len(connections_per_unique[side]) > 2: has_many = True
				for side2 in connections_per_unique[side]:
					if side2 in this_side: continue
					this_side.add(side2)
					added_any = True
					used_connections.add(canontip(side, side2))
					if len(connections_per_unique[side2]) > 2: has_many = True
		for side in this_side:
			for path in connections_per_unique[side]:
				if len(connections_per_unique[side]) == 2: has_two = True
				if len(connections_per_unique[side]) > 2: has_many = True
				for side2 in connections_per_unique[side]:
					if side2 in other_side: continue
					other_side.add(side2)
					added_any = True
					used_connections.add(canontip(side, side2))
					if len(connections_per_unique[side2]) > 2: has_many = True
		if not added_any: break
	if not has_two: return
	if has_many: return
	if len(other_side) != 2: return
	if len(this_side) != 2: return
	total_used_nodes = set()
	for node in this_side: total_used_nodes.add(node[1:])
	for node in other_side: total_used_nodes.add(node[1:])
	if len(total_used_nodes) != 4: return
	if len(used_connections) < 3: return
	if len(used_connections) > 4: return
	abundant_connections = list(used_connections)
	abundant_connections.sort(key=lambda x: -connection_counts[x])
	assert connection_counts[abundant_connections[1]] >= connection_counts[abundant_connections[2]]
	if connection_counts[abundant_connections[1]] == connection_counts[abundant_connections[2]]: return
	used_here = set()
	used_here.add(abundant_connections[0][0])
	used_here.add(abundant_connections[0][1])
	used_here.add(abundant_connections[1][0])
	used_here.add(abundant_connections[1][1])
	if len(used_here) < 4: return
	for i in range(2, len(abundant_connections)): forbidden_connections.add(abundant_connections[i])


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
for node in unique_nodes:
	check_side('>' + node)
	check_side('<' + node)

for pair in forbidden_connections:
	sys.stderr.write("removed connection: " + str(pair) + "\n")

for pathstr in all_paths:
	path = pathstr.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	if len(path) == 0: continue
	key = canontip(path[0], revnode(path[-1]))
	if key not in forbidden_connections:
		print(pathstr)
	else:
		sys.stderr.write("forbade " + pathstr + "\n")
