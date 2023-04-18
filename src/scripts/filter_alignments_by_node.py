#!/usr/bin/env python

import sys

read_file = sys.argv[1]
# alignments from stdin
# alignments to stdout

nodes= set()

with open(read_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		nodes.add(parts[0])

for l in sys.stdin:
	parts = l.strip().split('\t')
	aln_nodes=set(parts[5].replace(">", ",").replace("<", ",")[1:].split(","))
	
	if nodes.isdisjoint(aln_nodes) == False:
		print(l.strip()) 
