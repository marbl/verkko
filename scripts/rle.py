#!/usr/bin/env python

import fileinput

seq = ""
seq_name = ""
for l in fileinput.input():
	rle = l[0]
	if l[0] == '>':
		if len(seq) > 0:
			print(">" + seq_name)
			print(seq)
		seq = ""
		seq_name = l[1:].strip()
		continue
	l = l.strip().replace('a', 'A').replace('c', 'C').replace('g', 'G').replace('t', 'T')
	if len(seq) == 0:
		seq = l[0]
		l = l[1:]
	for c in l:
		if c != seq[-1]: seq += c

if len(seq) > 0:
	print(">" + seq_name)
	print(seq)
