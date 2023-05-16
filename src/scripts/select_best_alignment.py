#!/usr/bin/env python

import sys

alignment_paths = sys.argv[1]
winnowmap_alignment_file = sys.argv[2]
trimmed_alignment_file = sys.argv[3]
alignment_file = sys.argv[4]
min_mapq = int(sys.argv[5])
min_aln = int(sys.argv[6])

# alignments to stdout

# record paths from graph aligner
paths = {}
with open(alignment_paths) as f:
	for l in f:
		parts = l.strip().split('\t')
		n = "\t".join(set(parts[0].replace(">", " >").replace("<", " <").strip().split(' ')))
		#sys.stderr.write("Recording path %s for %s and current value is %s\n"%(parts[0], n, (paths[n] if n in paths else "")))
		paths[n] = parts[0]

# record where the read is mapped now
reads=set()
winnowmap_alns={}
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if (int(parts[3]) - int(parts[2]) < min_aln): continue
		if parts[0] not in winnowmap_alns:
			winnowmap_alns[parts[0]] = set()

		winnowmap_alns[parts[0]].add(parts[5])

# record where the read was mapped before (nodes)
alns={}
with open(trimmed_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] not in winnowmap_alns: continue
		if parts[0] not in alns:
			alns[parts[0]] = set()
		alns[parts[0]].update(parts[5].replace(">", " >").replace("<", " <").strip().split(' '))

# check every read we re-aligned, if it had no alignment before or the alignment changed, keep it
# otherwise keep the old version
todel = set()
for r in winnowmap_alns:
	if r not in alns: continue
	path = paths["\t".join(winnowmap_alns[r])] if "\t".join(winnowmap_alns[r]) in paths else set()

	#sys.stderr.write("Checking winnowmap aln for read %s which covers nodes %s and is equal %s and subset %s and path is %s\n"%(r, winnowmap_alns[r], alns[r] == winnowmap_alns[r], winnowmap_alns[r].issubset(alns[r]), path))
	# we keep a winnowmap read when it's not reproducing a known path and it either matches the same nodes as graph aligner (post triming) or it's it is not a subset of the graphaligner alignment nodes
	if path == set() and (alns[r] == winnowmap_alns[r] or winnowmap_alns[r].issubset(alns[r]) == False): continue
	todel.add(r)
[winnowmap_alns.pop(key) for key in todel]

# now output reads we didn't replace from old file
with open(alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		mapq=int(parts[11])
		if parts[0] not in winnowmap_alns and mapq >= min_mapq:
			print(l.strip())

# and new alignments
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		mapq=int(parts[11])
		if parts[0] in winnowmap_alns and mapq >= min_mapq:
			sys.stderr.write(parts[0] + "\n")
			print(l.strip())
