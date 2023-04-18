#!/usr/bin/env python

import sys

read_file = sys.argv[1]
min_mapq = int(sys.argv[2])
# alignments from stdin
# alignments to stdout

reads= set()

with open(read_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		reads.add(parts[0])

for l in sys.stdin:
	parts = l.strip().split('\t')
	mapq=int(parts[11])
	if parts[0] not in reads and mapq >= min_mapq:
		print(l.strip()) 
