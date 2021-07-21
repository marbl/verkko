#!/usr/bin/python

import sys

graph_file = sys.argv[1]
# paths gaf from stdin
# paths txt to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(e1, e2):
	fwstr = e1 + e2
	bwstr = e2 + e1
	if bwstr < fwstr: return (revnode(e2), revnode(e1))
	return (e1, e2)

node_exists = set()
edge_exists = set()
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S": node_exists.add(parts[1])
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			edge_exists.add(canon(fromnode, tonode))

for l in sys.stdin:
	parts = l.strip().split('\t')
	path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t')
	last_break = 0
	for i in range(0, len(path)):
		if path[i][1:] not in node_exists:
			if i > last_break: print("".join(path[last_break:i]))
			last_break = i+1
			continue
		if i > 0 and canon(path[i-1], path[i]) not in edge_exists:
			if i > last_break: print("".join(path[last_break:i]))
			last_break = i
			continue
	if last_break < len(path): print("".join(path[last_break:]))
