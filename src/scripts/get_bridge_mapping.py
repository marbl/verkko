#!/usr/bin/python

import sys

connected_graph_file = sys.argv[1]
not_connected_graph_file = sys.argv[2]
# mapping to stdout

existing_nodes = set()
with open(not_connected_graph_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			existing_nodes.add(parts[1])

with open(connected_graph_file) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			name = parts[1]
			if name not in existing_nodes:
				base_name = "_".join(name.split("_")[:-1])
				assert base_name in existing_nodes
				print(name + "\t" + ">" + base_name + ":0:0")
