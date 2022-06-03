#!/usr/bin/env python

import sys
import re

# gfa from stdin
# gfa to stdout

haplocut_prefix = "_hapcut"


# https://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
def atoi(text):
    return int(text) if text.isdigit() else text
def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]


nodes_per_haplocut = {}
existing_edges = set()
for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		name = parts[1]
		if haplocut_prefix in name:
			pos = name.find(haplocut_prefix)
			base_name = name[0:pos]
			if base_name not in nodes_per_haplocut: nodes_per_haplocut[base_name] = []
			nodes_per_haplocut[base_name].append(name)
	if parts[0] == "L":
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		existing_edges.add((fromnode, tonode))
	print(l.strip())

for basename in nodes_per_haplocut:
	assert len(nodes_per_haplocut[basename]) >= 2
	# the nodes might be further cut by chop_misassemblies.py
	# so we want to re-introduce edges in those cases as well
	# fortunately everything can be sorted by name as long as you take care of numbers
	nodes = nodes_per_haplocut[basename]
	nodes.sort(key=natural_keys)
	for i in range(1, len(nodes)):
		if (nodes[i-1], nodes[i]) in existing_edges: continue
		print("L\t" + nodes[i-1] + "\t+\t" + nodes[i] + "\t+\t0M")
