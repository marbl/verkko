#!/usr/bin/env python

import sys

relevant_reads_file = sys.argv[1]
old_gaf_file = sys.argv[2]
new_gaf_file = sys.argv[3]
# gaf to stdout

relevant_reads = set()
with open(relevant_reads_file) as f:
	for l in f:
		relevant_reads.add(l.strip())

with open(old_gaf_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0]
		if readname in relevant_reads: continue
		print(l.strip())

with open(new_gaf_file) as f:
	for l in f:
		print(l.strip())
