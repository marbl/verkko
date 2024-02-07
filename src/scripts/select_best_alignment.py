#!/usr/bin/env python

import sys
import statistics

def mad(values, mean):
	s = 0
	for i in values:
		s += abs(i-mean)
	return s/len(values)
 
alignment_paths = sys.argv[1]
winnowmap_alignment_file = sys.argv[2]
trimmed_alignment_file = sys.argv[3]
alignment_file = sys.argv[4]
min_mapq = int(sys.argv[5])
min_aln = int(sys.argv[6])
max_end_clip = 1500
idy_stdev = 1.5
graphaligner_winners = sys.argv[7]

graphaligner_tokeep = set()
with open(graphaligner_winners) as f:
	for l in f:
	    parts = l.strip().split('\t')
	    graphaligner_tokeep.add(parts[0])

# alignments to stdout

# record paths from graph aligner
paths = {}
with open(alignment_paths) as f:
	for l in f:
		parts = l.strip().split('\t')
		n = "\t".join(sorted(set(parts[0].replace(">", " >").replace("<", " <").strip().split(' '))))
		#sys.stderr.write("Recording path %s for %s and current value is %s\n"%(parts[0], n, (paths[n] if n in paths else "")))
		paths[n] = parts[0]

# record where the read is mapped now
currID=""
idys = []
read_alignment_idy = {}
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readend = int(parts[3])
		readstart = int(parts[2])
		leftclip = int(parts[7])
		rightclip = int(parts[6]) - int(parts[8])
		idy=float(parts[-1][5:])
		if (readend - readstart < min_aln): continue
		if (leftclip > max_end_clip and rightclip > max_end_clip): continue

		if  currID != parts[0]:
			if len(idys) > 0:
				assert(currID not in read_alignment_idy)
				read_alignment_idy[currID] = [statistics.median(idys), mad(idys, statistics.mean(idys))]
				#if len(idys) > 1: sys.stderr.write("Adding info for read %s which has mean %s and list %s and mad is %s sd %s\n"%(currID, statistics.median(idys), idys, mad(idys, statistics.mean(idys)), statistics.stdev(idys)))
				idys.clear()

		currID = parts[0]
		idys.append(idy)
	if len(idys) > 0:
		assert(currID not in read_alignment_idy)
		read_alignment_idy[currID] = [statistics.median(idys), mad(idys, statistics.mean(idys))]

# two pass, first records median alignment identity for a read, second only keeps those within some distand from median
reads=set()
winnowmap_alns={}
with open(winnowmap_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readend = int(parts[3])
		readstart = int(parts[2])
		leftclip = int(parts[7])
		rightclip = int(parts[6]) - int(parts[8])
		idy=float(parts[-1][5:])
		if (readend - readstart < min_aln): continue
		if (leftclip > max_end_clip and rightclip > max_end_clip): continue
		if read_alignment_idy[parts[0]][0] - idy_stdev * read_alignment_idy[parts[0]][1] > idy: continue

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
	path = paths["\t".join(sorted(winnowmap_alns[r]))] if "\t".join(sorted(winnowmap_alns[r])) in paths else set()

	#sys.stderr.write("Checking winnowmap aln for read %s which covers nodes %s vs %s and is equal %s and subset %s and path is %s\n"%(r, winnowmap_alns[r], alns[r], alns[r] == winnowmap_alns[r], winnowmap_alns[r].issubset(alns[r]), path))
	# we keep a winnowmap read when it's not reproducing a known path and it either matches the same nodes as graph aligner (post triming) or it's it is not a subset of the graphaligner alignment nodes
	if r not in graphaligner_tokeep and path == set() and (alns[r] == winnowmap_alns[r] or winnowmap_alns[r].issubset(alns[r]) == False): continue
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
		idy=float(parts[-1][5:])
		if parts[0] in winnowmap_alns:
			#sys.stderr.write("Looking up line %s and idy is %s vs expected %s\n"%(l.strip(), idy, (read_alignment_idy[parts[0]][0] - idy_stdev * read_alignment_idy[parts[0]][1])))
			if mapq >= min_mapq and idy >= read_alignment_idy[parts[0]][0] - idy_stdev * read_alignment_idy[parts[0]][1]:
				sys.stderr.write(parts[0] + "\n")
				print(l.strip())
