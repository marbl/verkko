#!/usr/bin/env python

import sys

uniquefile = sys.argv[1]
graphfile = sys.argv[2]
forbidden_connectionfile = sys.argv[3]
picked_connections_file = sys.argv[4]
pathsfile = sys.argv[5]
coveragefile = sys.argv[6]
min_solid_coverage = int(sys.argv[7])

def find(s, parent):
	while parent[s] != parent[parent[s]]:
		parent[s] = parent[parent[s]]
	return parent[s]

def union(s1, s2, parent, rank):
	p1 = find(s1, parent)
	p2 = find(s2, parent)
	if rank[p1] < rank[p2]: (p1, p2) = (p2, p1)
	parent[p2] = p1
	if rank[p2] == rank[p1]: rank[p1] += 1

def reverse(n):
	return (">" if n[0] == '<' else '<') + n[1:]

def canon(left, right):
	fwstr = left + right
	bwstr = reverse(right) + reverse(left)
	if bwstr < fwstr: return (reverse(right), reverse(left))
	return (left, right)

edge_coverage = {}
with open(pathsfile) as f:
	for l in f:
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		for i in range(1, len(path)):
			key = canon(path[i-1], path[i])
			if key not in edge_coverage: edge_coverage[key] = 0
			edge_coverage[key] += 1

solid_edges = set()
for key in edge_coverage:
	if edge_coverage[key] < min_solid_coverage: continue
	solid_edges.add(key)

unique_nodes = set()
with open(uniquefile) as f:
	for l in f:
		unique_nodes.add(l.strip())

solid_nodes = set()
with open(coveragefile) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[2] == "coverage": continue
		if float(parts[2]) < min_solid_coverage: continue
		if parts[0] in unique_nodes: continue
		solid_nodes.add(parts[0])

parent = {}
rank = {}

with open(graphfile) as f:
	for l in f:
		if l[0] == 'S':
			parts = l.strip().split('\t')
			if ">" + parts[1] not in parent: parent[">" + parts[1]] = ">" + parts[1]
			if "<" + parts[1] not in parent: parent["<" + parts[1]] = "<" + parts[1]
			if ">" + parts[1] not in rank: rank[">" + parts[1]] = 1
			if "<" + parts[1] not in rank: rank["<" + parts[1]] = 1
			if parts[1] not in unique_nodes:
				if parts[1] not in parent: parent[parts[1]] = parts[1]
				if parts[1] not in rank: rank[parts[1]] = 1
				union(">" + parts[1], parts[1], parent, rank)
				union(">" + parts[1], "<" + parts[1], parent, rank)
		if l[0] == 'L':
			parts = l.strip().split('\t')
			fromend = (">" if parts[2] == '+' else '<') + parts[1]
			toend = ("<" if parts[4] == '+' else '>') + parts[3]
			if fromend not in parent: parent[fromend] = fromend
			if toend not in parent: parent[toend] = toend
			if fromend not in rank: rank[fromend] = 1
			if toend not in rank: rank[toend] = 1
			union(fromend, toend, parent, rank)

edge_needs_connection = set(solid_edges)
node_needs_connection = set(solid_nodes)
has_connection = set()

sys.stderr.write(str(len(node_needs_connection)) + " nodes need covering\n")
sys.stderr.write(str(len(edge_needs_connection)) + " edges need covering\n")

with open(forbidden_connectionfile) as f:
	for line in f:
		l = line.strip()
		last_break = 0
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		assert len(path) >= 2
		for node in path:
			if node[1:] in node_needs_connection: node_needs_connection.remove(node[1:])
		for i in range(1, len(path)):
			key = canon(path[i-1], path[i])
			if key in edge_needs_connection: edge_needs_connection.remove(key)
		fwkey = path[0]
		bwkey = reverse(path[-1])
		has_connection.add(fwkey)
		has_connection.add(bwkey)

with open(picked_connections_file) as f:
	for line in f:
		l = line.strip()
		last_break = 0
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		assert len(path) >= 2
		for node in path:
			if node[1:] in node_needs_connection: node_needs_connection.remove(node[1:])
		for i in range(1, len(path)):
			key = canon(path[i-1], path[i])
			if key in edge_needs_connection: edge_needs_connection.remove(key)
		fwkey = path[0]
		bwkey = reverse(path[-1])
		has_connection.add(fwkey)
		has_connection.add(bwkey)

sys.stderr.write(str(len(node_needs_connection)) + " uncovered nodes\n")
sys.stderr.write(str(len(edge_needs_connection)) + " uncovered edges\n")

forbidden_tangles = set()

for node in node_needs_connection:
	if node not in parent: continue
	forbidden_tangles.add(find(node, parent))
for edge in edge_needs_connection:
	if edge[0] not in parent or edge[1] not in parent: continue
	forbidden_tangles.add(find(edge[0], parent))

for node in unique_nodes:
	if ">" + node not in has_connection:
		forbidden_tangles.add(find(">" + node, parent))
	if "<" + node not in has_connection:
		forbidden_tangles.add(find("<" + node, parent))

for node in unique_nodes:
	if find(">" + node, parent) in forbidden_tangles:
		print(">" + node)
	if find("<" + node, parent) in forbidden_tangles:
		print("<" + node)
