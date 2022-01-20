#!/usr/bin/env python

import sys

base_graph = sys.argv[1]
unrolled_graph = sys.argv[2]
# mapping to stdout

existing_nodes = set()
with open(base_graph) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			existing_nodes.add(parts[1])

with open(unrolled_graph) as f:
	for l in f:
		parts = l.strip().split("\t")
		if parts[0] == "S":
			name = parts[1]
			if name not in existing_nodes:
				assert name[:7] == "unroll_"
				base_name = "_".join(name[7:].split("_")[:-1])
				assert base_name in existing_nodes
				print(name + "\t" + ">" + base_name + ":0:0")
