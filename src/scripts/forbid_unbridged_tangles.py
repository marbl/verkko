#!/usr/bin/env python

import sys
import graph_functions as gf

uniquefile = sys.argv[1]
graphfile = sys.argv[2]
forbidden_connectionfile = sys.argv[3]
picked_connections_file = sys.argv[4]
pathsfile = sys.argv[5]
ontcoveragefile = sys.argv[6]
min_ont_solid_coverage = int(sys.argv[7])
hificoveragefile = sys.argv[8]
min_hifi_solid_coverage = int(sys.argv[9])


length_solid_node_threshold = 200000

def union(s1, s2, parent, rank):
	p1 = gf.find(parent, s1)
	p2 = gf.find(parent, s2)
	if rank[p1] < rank[p2]: (p1, p2) = (p2, p1)
	parent[p2] = p1
	if rank[p2] == rank[p1]: rank[p1] += 1

edge_coverage = {}
with open(pathsfile) as f:
	for l in f:
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		for i in range(1, len(path)):
			key = gf.canon(path[i-1], path[i])
			if key not in edge_coverage: edge_coverage[key] = 0
			edge_coverage[key] += 1

unique_nodes = set()
with open(uniquefile) as f:
	for l in f:
		unique_nodes.add(l.strip())

hifi_solid_nodes = set()
with open(hificoveragefile) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node":
			assert(parts[1] == "coverage" and parts[2] == "length")
			continue
		if float(parts[1]) < min_hifi_solid_coverage: continue
		if parts[0] in unique_nodes: continue
		hifi_solid_nodes.add(parts[0])

length_solid_nodes = set()
ont_solid_nodes = set()
with open(ontcoveragefile) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node":
			assert(parts[1] == "coverage" and parts[2] == "length")
			continue
		if parts[0] in unique_nodes: continue
		if float(parts[1]) >= min_ont_solid_coverage: ont_solid_nodes.add(parts[0])
		if int(parts[2]) >= length_solid_node_threshold: length_solid_nodes.add(parts[0])

solid_nodes = set()
for node in hifi_solid_nodes:
	if node in ont_solid_nodes:
		solid_nodes.add(node)
for node in length_solid_nodes:
	solid_nodes.add(node)

solid_edges = set()
for key in edge_coverage:
	if edge_coverage[key] < min_ont_solid_coverage: continue
	if key[0][1:] not in hifi_solid_nodes or key[1][1:] not in hifi_solid_nodes: continue
	solid_edges.add(key)

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
			key = gf.canon(path[i-1], path[i])
			if key in edge_needs_connection: edge_needs_connection.remove(key)
		fwkey = path[0]
		bwkey = gf.revnode(path[-1])
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
			key = gf.canon(path[i-1], path[i])
			if key in edge_needs_connection: edge_needs_connection.remove(key)
		fwkey = path[0]
		bwkey = gf.revnode(path[-1])
		has_connection.add(fwkey)
		has_connection.add(bwkey)

sys.stderr.write(str(len(node_needs_connection)) + " uncovered nodes\n")
sys.stderr.write(str(len(edge_needs_connection)) + " uncovered edges\n")

forbidden_tangles = set()

for node in node_needs_connection:
	if node not in parent: continue
	forbidden_tangles.add(gf.find(parent, node))
for edge in edge_needs_connection:
	if edge[0] not in parent or edge[1] not in parent: continue
	forbidden_tangles.add(gf.find(parent, edge[0]))

for node in unique_nodes:
	if ">" + node not in has_connection:
		forbidden_tangles.add(gf.find(parent, ">" + node))
	if "<" + node not in has_connection:
		forbidden_tangles.add(gf.find(parent, "<" + node))

for node in unique_nodes:
	if gf.find(parent, ">" + node) in forbidden_tangles:
		print(">" + node)
	if gf.find(parent, "<" + node) in forbidden_tangles:
		print("<" + node)
