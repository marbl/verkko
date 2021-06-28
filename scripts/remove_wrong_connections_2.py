#!/usr/bin/python

import fileinput

def pathstr(p):
	return "".join(p)

def reverse(n):
	return (">" if n[0] == '<' else '<') + n[1:]

connectors = {}

connections = []

for line in fileinput.input():
	l = line.strip() + ">"
	last_break = 0
	path = []
	for i in range(1, len(l)):
		if l[i] == '<' or l[i] == '>':
			path.append(l[last_break:i])
			last_break = i
	assert len(path) >= 2
	fwkey = path[0]
	bwkey = reverse(path[-1])
	if fwkey not in connectors: connectors[fwkey] = {}
	if bwkey not in connectors[fwkey]: connectors[fwkey][bwkey] = []
	connectors[fwkey][bwkey].append(len(connections))
	if bwkey not in connectors: connectors[bwkey] = {}
	if fwkey not in connectors[bwkey]: connectors[bwkey][fwkey] = []
	connectors[bwkey][fwkey].append(len(connections))
	connections.append(path)

forbidden = set()

for fwpos in connectors:
	max_coverage = 0
	for bwpos in connectors[fwpos]:
		max_coverage = max(max_coverage, len(connectors[fwpos][bwpos]))
	for bwpos in connectors[fwpos]:
		if max_coverage >= 2 * len(connectors[fwpos][bwpos]):
			for j in connectors[fwpos][bwpos]:
				forbidden.add(j)

for i in range(0, len(connections)):
	if i in forbidden: continue
	print(pathstr(connections[i]))
