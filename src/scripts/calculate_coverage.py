#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
# gaf from stdin

node_sizes = {}

with open(graph_file) as f:
	for l in f:
		if l[0] != 'S': continue
		parts = l.strip().split('\t')
		node_sizes[parts[1]] = len(parts[2])

node_average_coverages = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if len(parts) < 10: continue
	pathstr = parts[5] + '>'
	pathsize = int(parts[6])
	pathstart = int(parts[7])
	pathend = int(parts[8])
	left_clip = pathstart
	right_clip = pathsize - pathend
	last_break = 0
	path = []
	for i in range(1, len(pathstr)):
		if pathstr[i] == '<' or pathstr[i] == '>':
			path.append(pathstr[last_break:i])
			last_break = i
	if len(path) == 1:
		nodename = path[0][1:]
		assert nodename in node_sizes
		if nodename not in node_average_coverages: node_average_coverages[nodename] = 0.0
		node_average_coverages[nodename] += float(pathsize - left_clip - right_clip) / float(pathsize)
		continue
	for i in range(1, len(path)-1):
		nodename = path[i][1:]
		# skip gap nodes that aren't in the graph and have no coverage to compute anyway
		if nodename not in node_sizes:
			continue
		if nodename not in node_average_coverages: node_average_coverages[nodename] = 0.0
		node_average_coverages[nodename] += 1.0
	# unlike above, we do want to make sure the path starts and ends with a valid node or else we have a problem
	nodename = path[0][1:]
	assert nodename in node_sizes
	if nodename not in node_average_coverages: node_average_coverages[nodename] = 0.0
	node_average_coverages[nodename] += float(node_sizes[nodename] - left_clip) / float(node_sizes[nodename])
	nodename = path[-1][1:]
	assert nodename in node_sizes
	if nodename not in node_average_coverages: node_average_coverages[nodename] = 0.0
	node_average_coverages[nodename] += float(node_sizes[nodename] - right_clip) / float(node_sizes[nodename])

print("node\tcoverage\tlength")
for node in node_sizes:
	cov = 0.0
	if node in node_average_coverages: cov = node_average_coverages[node]
	print(node + "\t" + str(cov) + "\t" + str(node_sizes[node]))
