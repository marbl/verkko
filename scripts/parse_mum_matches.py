#!/usr/bin/python

import sys
import re

kmer_size = int(sys.argv[1])
# mums from stdin

match_count = {}
current_contig = ""
for l in sys.stdin:
	if l[0] == '>':
		new_contig = l[1:].strip()
		if len(new_contig) > 8 and new_contig[-8:] == " Reverse": new_contig = new_contig[:-8]
		if current_contig != new_contig:
			print(current_contig + "\t" + "\t".join(str(key) + ":" + str(match_count[key]) for key in match_count))
			match_count = {}
		current_contig = new_contig
		continue
	parts = re.sub("\\s+", "\t", l.strip()).split('\t')
	match_len = int(parts[3])
	name = parts[0][0:5]
	if match_len < kmer_size: continue
	if name not in match_count: match_count[name] = 0
	match_count[name] += match_len - kmer_size + 1

print(current_contig + "\t" + "\t".join(str(key) + ":" + str(match_count[key]) for key in match_count) + "\t" + ":".join(str(match_count[key]) for key in match_count))
