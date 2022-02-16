#!/usr/bin/env python

import sys

prefix = sys.argv[1]
mapping_file = sys.argv[2]
# gfa from stdin
# gfa to stdout

def iterate_deterministic(l):
	tmp = list(l)
	tmp.sort()
	for x in tmp:
		yield x

def getone(s):
	assert len(s) == 1
	for n in s:
		return n

# https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
trans = str.maketrans("ACTG", "TGAC")
def revcomp(s):
	return s.translate(trans)[::-1]

def revnode(n):
	return (">" if n[0] == "<" else "<") + n[1:]

def start_unitig(startpos, unitigs, belongs_to_unitig, edges):
	assert startpos[1:] not in belongs_to_unitig
	new_unitig = [startpos]
	belongs_to_unitig.add(startpos[1:])
	while len(edges[new_unitig[-1]]) == 1 and getone(edges[new_unitig[-1]])[1:] != new_unitig[-1][1:]:
		new_pos = getone(edges[new_unitig[-1]])
		if len(edges[revnode(new_pos)]) != 1: break
		assert new_pos[1:] not in belongs_to_unitig
		new_unitig.append(new_pos)
		belongs_to_unitig.add(new_pos[1:])
	unitigs.append(new_unitig)

def start_circular_unitig(startpos, unitigs, belongs_to_unitig, edges):
	assert startpos[1:] not in belongs_to_unitig
	assert len(edges[revnode(startpos)]) == 1
	assert len(edges[startpos]) == 1
	new_unitig = [startpos]
	belongs_to_unitig.add(startpos[1:])
	while True:
		assert len(edges[new_unitig[-1]]) == 1
		new_pos = getone(edges[new_unitig[-1]])
		assert len(edges[revnode(new_pos)]) == 1
		if new_pos == startpos: break
		assert new_pos[1:] not in belongs_to_unitig
		new_unitig.append(new_pos)
		belongs_to_unitig.add(new_pos[1:])
	unitigs.append(new_unitig)

def write_seq(stream, unitig, node_seqs, edge_overlaps):
	seq = node_seqs[unitig[0][1:]]
	if unitig[0][0] == "<": seq = revcomp(seq)
	stream.write(seq)
	for i in range(1, len(unitig)):
		add = node_seqs[unitig[i][1:]]
		if unitig[i][0] == "<": add = revcomp(add)
		stream.write(add[edge_overlaps[(unitig[i-1], unitig[i])]:])


node_seqs = {}
edges = {}
edge_overlaps = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		node_seqs[parts[1]] = parts[2]
		if ">" + parts[1] not in edges: edges[">" + parts[1]] = set()
		if "<" + parts[1] not in edges: edges["<" + parts[1]] = set()
	elif parts[0] == 'L':
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		edges[fromnode].add(tonode)
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[revnode(tonode)].add(revnode(fromnode))
		edge_overlaps[(fromnode, tonode)] = int(parts[5][:-1])
		edge_overlaps[(revnode(tonode), revnode(fromnode))] = int(parts[5][:-1])

belongs_to_unitig = set()
unitigs = []

nodenames = list(node_seqs)
nodenames.sort()

for node in nodenames:
	assert ">" + node in edges
	assert "<" + node in edges
	if len(edges[">" + node]) != 1:
		if node not in belongs_to_unitig: start_unitig("<" + node, unitigs, belongs_to_unitig, edges)
		for edge in iterate_deterministic(edges[">" + node]):
			if edge[1:] not in belongs_to_unitig: start_unitig(edge, unitigs, belongs_to_unitig, edges)
	if len(edges["<" + node]) != 1:
		if node not in belongs_to_unitig: start_unitig(">" + node, unitigs, belongs_to_unitig, edges)
		for edge in iterate_deterministic(edges["<" + node]):
			if edge[1:] not in belongs_to_unitig: start_unitig(edge, unitigs, belongs_to_unitig, edges)
	if len(edges[">" + node]) == 1 and getone(edges[">" + node])[1:] == node:
		if node not in belongs_to_unitig: start_unitig("<" + node, unitigs, belongs_to_unitig, edges)
	if len(edges["<" + node]) == 1 and getone(edges["<" + node])[1:] == node:
		if node not in belongs_to_unitig: start_unitig(">" + node, unitigs, belongs_to_unitig, edges)

for node in nodenames:
	if node in belongs_to_unitig: continue
	start_circular_unitig(">" + node, unitigs, belongs_to_unitig, edges)

unitig_start = {}
unitig_end = {}
for i in range(0, len(unitigs)):
	unitig_start[unitigs[i][0]] = ">" + str(i)
	unitig_start[revnode(unitigs[i][-1])] = "<" + str(i)
	unitig_end[unitigs[i][-1]] = ">" + str(i)
	unitig_end[revnode(unitigs[i][0])] = "<" + str(i)

unitig_edges = set()
for node in nodenames:
	assert node in belongs_to_unitig
	if ">" + node in unitig_end:
		for target in iterate_deterministic(edges[">" + node]):
			assert target in unitig_start
			edge = (unitig_end[">" + node], unitig_start[target], edge_overlaps[(">" + node, target)])
			unitig_edges.add(edge)
	if "<" + node in unitig_end:
		for target in iterate_deterministic(edges["<" + node]):
			assert target in unitig_start
			edge = (unitig_end["<" + node], unitig_start[target], edge_overlaps[("<" + node, target)])
			unitig_edges.add(edge)

with open(mapping_file, "w") as f:
	for i in range(0, len(unitigs)):
		f.write(prefix + str(i) + "\t" + "".join(unitigs[i]) + ":0:0" + "\n")

for i in range(0, len(unitigs)):
	sys.stdout.write("S\t" + prefix + str(i) + "\t")
	write_seq(sys.stdout, unitigs[i], node_seqs, edge_overlaps)
	sys.stdout.write("\n")

unitig_edges = list(unitig_edges)
unitig_edges.sort()

for edge in unitig_edges:
	print("L\t" + prefix + edge[0][1:] + "\t" + ("+" if edge[0][0] == ">" else "-") + "\t" + prefix + edge[1][1:] + "\t" + ("+" if edge[1][0] == ">" else "-") + "\t" + str(edge[2]) + "M")
