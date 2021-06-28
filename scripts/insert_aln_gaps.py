#!/usr/bin/python

import sys

graph_file = sys.argv[1]
min_gap_coverage = int(sys.argv[2])
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
alns = []
last_read_name = ""

for l in sys.stdin:
	parts = l.strip().split('\t')
	readname = parts[0].split(' ')[0]
	if readname != last_read_name:
		if len(alns) >= 2:
			alns.sort(key=lambda x: x[0])
			for i in range(1, len(alns)):
				assert alns[i][0] >= alns[i-1][0]
				gap_len = alns[i][0] - alns[i-1][1]
				if alns[i-1][3] in not_tips or alns[i][2] in not_tips: continue
				if alns[i-1][5] > 50 or alns[i][4] > 50: continue
				key = canontip(alns[i-1][3], alns[i][2])
				if key not in gaps: gaps[key] = []
				gaps[key].append(gap_len)
		alns = []
		last_read_name = readname
	readstart = int(parts[2])
	readend = int(parts[3])
	path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t')
	alnstart = revnode(path[0])
	alnend = path[-1]
	leftclip = int(parts[7])
	rightclip = int(parts[6]) - int(parts[8])
	alns.append((readstart, readend, alnstart, alnend, leftclip, rightclip))

next_gap_id = 1

for gap in gaps:
	if len(gaps[gap]) < min_gap_coverage: continue
	gap_len_sum = 0
	gap_len_count = 0
	for length in gaps[gap]:
		gap_len_sum += length
		gap_len_count += 1
	gap_len = gap_len_sum / gap_len_count
	gap_name = "gap-" + str(next_gap_id) + "-len-" + str(gap_len) + "-cov-" + str(len(gaps[gap]))
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
	if gap_len < 1: gap_len = 1
	next_gap_id += 1
	overlap = str(overlap) + "M"
	print("S\t" + gap_name + "\t" + seq)
	print("L\t" + gap[0][1:] + "\t" + ("+" if gap[0][0] == ">" else "-") + "\t" + gap_name + "\t+\t" + overlap)
	print("L\t" + gap_name + "\t+\t" + gap[1][1:] + "\t" + ("-" if gap[1][0] == ">" else "+") + "\t" + overlap)

sys.stderr.write("inserted " + str(next_gap_id-1) + " gaps\n")
