#!/usr/bin/env python

import sys

winnowmap_alignment_file = sys.argv[1]
alignment_file = sys.argv[2]
min_mapq = int(sys.argv[3])
# alignments to stdout

# record where the read is mapped now
reads=set()
winnowmap_alns={}
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] not in winnowmap_alns:
			winnowmap_alns[parts[0]] = set()
		winnowmap_alns[parts[0]].add(parts[5])
		sys.stderr.write("Winnowmap adding read %s with part %s\n"%(parts[0], winnowmap_alns[parts[0]]))

# record where the read was mapped before (nodes)
alns={}
with open(alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] not in winnowmap_alns: continue

		if parts[0] not in alns:
			alns[parts[0]] = set()
		alns[parts[0]].update(parts[5].replace(">", " >").replace("<", " <").strip().split(' '))
		sys.stderr.write("Adding read %s with nodes %s\n"%(parts[0], alns[parts[0]]))

# check every read we re-aligned, if it had no alignment before or the alignment changed, keep it
# otherwise keep the old version
todel = set()
for r in winnowmap_alns:
	if r not in alns:
		continue
	sys.stderr.write("For read %s comparing winnowmap %s to graphaligner %s shared is %s subset is %s equal is %s\n"%(r, winnowmap_alns[r], alns[r], winnowmap_alns[r].intersection(alns[r]),winnowmap_alns[r].issubset(alns[r]),winnowmap_alns[r] == alns[r] ))
	if alns[r] == winnowmap_alns[r] or winnowmap_alns[r].issubset(alns[r]) == False:
		sys.stderr.write("Here keeping read %s because the if is %s\n"%(r, alns[r] == winnowmap_alns[r] or winnowmap_alns[r].issubset(alns[r]) == False))
		continue
	sys.stderr.write("Removing read %s\n"%(r))
	todel.add(r)
[winnowmap_alns.pop(key) for key in todel]

# now output reads we didn't replace from old file
with open(alignment_file) as f:
	for l in f:	
		parts = l.strip().split('\t')
		mapq=int(parts[11])
		if parts[0] not in winnowmap_alns and mapq >= min_mapq:
			print(l.strip()) 
			sys.stderr.write("Writing graph aligner for read %s because it is not in the set to output\n"%(parts[0]))

# and new alignments
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		mapq=int(parts[11])
		if parts[0] in winnowmap_alns and mapq >= min_mapq:
			sys.stderr.write("Writing winnowmap alignment for read %s because it is in the ste to output\n"%(parts[0]))
			print(l.strip())
