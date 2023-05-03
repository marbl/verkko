#!/usr/bin/env python

import sys

def str2bool(v):
  return v.lower() in ("yes", "true", "t", "1")


read_file = sys.argv[1]
filter_offset = int(sys.argv[2])
invert = str2bool(sys.argv[3])
# alignments from stdin
# alignments to stdout

nodes= set()

with open(read_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		nodes.add(parts[0])

for l in sys.stdin:
	parts = l.strip().split('\t')
	if ">" in parts[filter_offset]:
		aln_nodes=set(parts[filter_offset].replace(">", ",").replace("<", ",")[1:].split(","))
	else:
		aln_nodes=set(parts[filter_offset].split(","))
	
	if invert == False and nodes.isdisjoint(aln_nodes) == False:
		print(l.strip()) 
	elif invert == True and nodes.isdisjoint(aln_nodes) == True:
		print(l.strip())
