#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
min_gap_coverage = int(sys.argv[2])
max_end_clip = int(sys.argv[3])
gap_aln_file_out = sys.argv[4]
prefix = sys.argv[5]
# gaf from stdin
# graph to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canontip(left, right):
	fwstr = left + right
	bwstr = right + left
	if bwstr < fwstr: return (right, left)
	return (left, right)

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
	return "".join(comp[c] for c in s[::-1])

not_tips = set()
node_seqs = {}

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			node_seqs[parts[1]] = parts[2]
		if parts[0] == 'L':
			not_tips.add((">" if parts[2] == "+" else "<") + parts[1])
			not_tips.add(("<" if parts[4] == "+" else ">") + parts[3])
		print(l.strip())

gaps = {}
reads_per_gap = {}
alns = []
last_read_name = ""

alns_per_read = {}
for l in sys.stdin:
	parts = l.strip().split('\t')
	readname = parts[0].split(' ')[0]
	readlen = int(parts[1])
	readstart = int(parts[2])
	readend = int(parts[3])
	path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t')
	alnstart = revnode(path[0])
	alnend = path[-1]
	leftclip = int(parts[7])
	rightclip = int(parts[6]) - int(parts[8])
	if readname not in alns_per_read: alns_per_read[readname] = []
	alns_per_read[readname].append((readstart, readend, alnstart, alnend, leftclip, rightclip, readlen))

for name in alns_per_read:
	alns = alns_per_read[name]
	if len(alns) >= 2:
		alns.sort(key=lambda x: x[0])
		for i in range(1, len(alns)):
			assert alns[i][0] >= alns[i-1][0]
			gap_len = (alns[i][0] - alns[i][4]) - (alns[i-1][1] + alns[i-1][5])
			if alns[i-1][3] in not_tips or alns[i][2] in not_tips: continue
			if alns[i-1][5] > max_end_clip or alns[i][4] > max_end_clip: continue
			key = canontip(alns[i-1][3], alns[i][2])
			fw = (key == (alns[i-1][3], alns[i][2]))
			if key not in gaps: gaps[key] = []
			gaps[key].append(gap_len)
			if key not in reads_per_gap: reads_per_gap[key] = set()
			readlen = alns[i][6]
			readstart = alns[i-1][1] + alns[i-1][5]
			readend = alns[i][0] - alns[i][4]
			reads_per_gap[key].add((name, readlen, readstart, readend, fw))

next_gap_id = 1
used_reads = set()
gap_lens = {}

for gap in gaps:
	if len(gaps[gap]) < min_gap_coverage: continue
	gap_len_list = list(gaps[gap])
	gap_len_list.sort()
	gap_len = gap_len_list[len(gap_len_list)/2]
	gap_name = prefix + "-" + str(next_gap_id) + "-len-" + str(gap_len) + "-cov-" + str(len(gaps[gap]))
	overlap = 0
	seq = "N"
	if gap_len < 0:
		overlap = -gap_len
		if overlap >= len(node_seqs[gap[0][1:]])-1: overlap = len(node_seqs[gap[0][1:]])-2
		if overlap >= len(node_seqs[gap[1][1:]])-1: overlap = len(node_seqs[gap[1][1:]])-2
		base_seq = node_seqs[gap[0][1:]]
		if gap[0][0] == "<": base_seq = revcomp(base_seq)
		assert len(base_seq) > overlap+1
		seq = base_seq[-overlap-1:]
		base_seq = node_seqs[gap[1][1:]]
		if gap[1][0] == ">": base_seq = revcomp(base_seq)
		assert len(base_seq) > overlap+1
		seq += base_seq[overlap]
	elif gap_len == 0:
		base_seq = node_seqs[gap[0][1:]]
		if gap[0][0] == "<": base_seq = revcomp(base_seq)
		seq = base_seq[-1]
		base_seq = node_seqs[gap[1][1:]]
		if gap[1][0] == ">": base_seq = revcomp(base_seq)
		seq += base_seq[0]
		overlap = 1
	else:
		seq = "N" * gap_len
		overlap = 0
	gap_lens[gap_name] = len(seq)
	if gap_len < 1: gap_len = 1
	next_gap_id += 1
	overlap = str(overlap) + "M"
	print("S\t" + gap_name + "\t" + seq)
	print("L\t" + gap[0][1:] + "\t" + ("+" if gap[0][0] == ">" else "-") + "\t" + gap_name + "\t+\t" + overlap)
	print("L\t" + gap_name + "\t+\t" + gap[1][1:] + "\t" + ("-" if gap[1][0] == ">" else "+") + "\t" + overlap)
	for name in reads_per_gap[gap]: used_reads.add((name, gap_name))

with open(gap_aln_file_out, "w") as f:
	for pair in used_reads:
		name = pair[0]
		gap_name = pair[1]
		f.write(name[0] + "\t" + str(name[1]) + "\t" + str(name[2]) + "\t" + str(name[3]) + "\t+\t" + (">" if name[4] else "<") + gap_name + "\t" + str(gap_lens[gap_name]) + "\t0\t" + str(gap_lens[gap_name]) + "\t" + str(gap_lens[gap_name]) + "\t" + str(gap_lens[gap_name]) + "\t" + "60" + "\n")

sys.stderr.write("inserted " + str(next_gap_id-1) + " gaps\n")
