#!/usr/bin/env python

import sys
import graph_functions as gf

forbidden_paths_out_file = sys.argv[1]
# paths from stdin
# allowed paths to stdout

def canon(p):
	ps = gf.pathstr(p)
	reverse = [('>' if n[0] == '<' else '<') + n[1:] for n in p[::-1]]
	rps = gf.pathstr(reverse)
	if ps < rps: return p
	return reverse

connectors = {}

for line in sys.stdin:
	l = line.strip() + ">"
	last_break = 0
	path = []
	for i in range(1, len(l)):
		if l[i] == '<' or l[i] == '>':
			path.append(l[last_break:i])
			last_break = i
	assert len(path) >= 2
	key = tuple(canon([path[0], path[-1]]))
	path = tuple(canon(path))
	if key not in connectors: connectors[key] = {}
	if path not in connectors[key]: connectors[key][path] = 0
	connectors[key][path] += 1

forbidden = set()

for key in connectors:
	connectorlist = []
	for path in connectors[key]:
		connectorlist.append((path, connectors[key][path]))
	connectorlist.sort(key=lambda x: -x[1])
	assert len(connectorlist) >= 1
	for connection in connectorlist:
		if connection[1]*2 >= connectorlist[0][1]:
			print(gf.pathstr(list(connection[0])))
		else:
			forbidden.add(gf.pathstr(list(connection[0])))

with open(forbidden_paths_out_file, "w") as f:
	for path in forbidden:
		f.write(path + "\n")

