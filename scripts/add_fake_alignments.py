#!/usr/bin/python

import sys

input_raw_graph_file = sys.argv[1]
input_connected_graph_file = sys.argv[2]
input_alignments_file = sys.argv[3]
input_node_coverages_csv_file = sys.argv[4]
output_fake_alignments_file = sys.argv[5]
output_fake_node_coverages_file = sys.argv[6]
fake_coverage = int(sys.argv[7])

real_nodes = set()
with open(input_raw_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			real_nodes.add(parts[1])

fake_nodes = {}
fake_edges = set()
with open(input_connected_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if parts[1] not in real_nodes:
				fake_nodes[parts[1]] = len(parts[2])
		if parts[0] == 'L':
			if parts[1] not in real_nodes or parts[3] not in real_nodes:
				fake_edges.add(((">" if parts[2] == "+" else "<") + parts[1], (">" if parts[4] == "+" else "<") + parts[3]))


with open(output_fake_node_coverages_file, "w") as f:
	with open(input_node_coverages_csv_file) as f2:
		for l in f2:
			f.write(l)
	for node in fake_nodes:
		f.write(node + "\t" + str(fake_nodes[node]) + "\t" + str(int(fake_nodes[node]) * fake_coverage) + "\n")

with open(output_fake_alignments_file, "w") as f:
	with open(input_alignments_file) as f2:
		for l in f2:
			f.write(l)
	for edge in fake_edges:
		for i in range(0, fake_coverage):
			f.write("FAKE\t1\t0\t1\t+\t" + edge[0] + edge[1] + "\t1\t0\t1\t1\t1\t60\tNM:i:0\tAS:f:100000\tdv:f:0\tid:f:1\n")


