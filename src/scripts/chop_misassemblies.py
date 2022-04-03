#!/usr/bin/env python

import sys

in_graph_file = sys.argv[1]
in_alns_file = sys.argv[2]
out_mapping_file = sys.argv[3]
out_cut_alns_file = sys.argv[4]
cut_coverage_threshold = int(sys.argv[5])
cut_length_threshold = int(sys.argv[6])
# graph to stdout

cut_anchor_min_size = 5000
min_dist_from_end = 50

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:]

nodelens = {}
max_overlap = {}
with open(in_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodelens[parts[1]] = len(parts[2])
		elif parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			overlap = int(parts[5][:-1])
			if fromnode not in max_overlap: max_overlap[fromnode] = 0
			if revnode(tonode) not in max_overlap: max_overlap[revnode(tonode)] = 0
			max_overlap[fromnode] = max(max_overlap[fromnode], overlap)
			max_overlap[revnode(tonode)] = max(max_overlap[revnode(tonode)], overlap)

alns_per_read = {}
with open(in_alns_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0]
		if readname not in alns_per_read: aln_end_positions[readname] = 0
		alns_per_read[readname] += 1

read_aln_positions = {}
with open(in_alns_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0]
		assert alns_per_read[readname] >= 1
		if alns_per_read[readname] == 1: continue
		path = parts[5].replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		alnstart = int(parts[2])
		alnend = int(parts[3])
		start_node_offset = int(parts[7])
		start_node = path[0][1:]
		if path[0][0] == "<": start_node_offset = nodelens[start_node] - start_node_offset - 1
		end_node = path[-1][1:]
		end_node_offset = int(parts[6]) - int(parts[8])
		if path[-1][0] == ">": end_node_offset = nodelens[end_node] - end_node_offset - 1
		if readname not in read_aln_positions: read_aln_positions[readname] = []
		read_aln_positions[readname].append((alnstart, alnend, start_node, start_node_offset, end_node, end_node_offset))

aln_end_positions = {}
for read in read_aln_positions:
	assert len(read_aln_positions[read]) >= 2
	read_aln_positions[read].sort(key= lambda x: x[0])
	for i in range(1, len(read_aln_positions[read])):
		assert read_aln_positions[read][i][0] >= read_aln_positions[read][i-1][0]
		if read_aln_positions[read][i][0] == read_aln_positions[read][i-1][0]: continue
		if read_aln_positions[read][i][1] <= read_aln_positions[read][i-1][1]: continue
		prev_node = read_aln_positions[read][i-1][4]
		this_node = read_aln_positions[read][i][2]
		prev_node_offset = read_aln_positions[read][i-1][5]
		this_node_offset = read_aln_positions[read][i][3]
		if prev_node not in aln_end_positions: aln_end_positions[prev_node] = []
		if this_node not in aln_end_positions: aln_end_positions[this_node] = []
		aln_end_positions[prev_node].append(prev_node_offset)
		aln_end_positions[this_node].append(this_node_offset)

cut_positions = {}
for node in aln_end_positions:
	aln_end_positions[node].sort()
	cut_start = None
	for i in range(0, len(aln_end_positions[node]) - cut_coverage_threshold):
		if aln_end_positions[node][i+cut_coverage_threshold] - aln_end_positions[node][i] <= cut_length_threshold:
			if cut_start is None:
				cut_start = i
			continue
		if cut_start is not None:
			start = cut_start
			end = i - 1 + cut_coverage_threshold
			mid = int((end - start) / 2) + start
			if node not in cut_positions: cut_positions[node] = []
			if len(cut_positions[node]) == 0 or cut_positions[node][-1] != aln_end_positions[node][mid]:
				cut_positions[node].append(aln_end_positions[node][mid])
		cut_start = None
	if cut_start is not None:
		start = cut_start
		end = len(aln_end_positions[node])-1
		mid = int((end - start) / 2) + start
		if node not in cut_positions: cut_positions[node] = []
		if len(cut_positions[node]) == 0 or cut_positions[node][-1] != aln_end_positions[node][mid]:
			cut_positions[node].append(aln_end_positions[node][mid])

for node in cut_positions:
	if len(cut_positions[node]) < 1: continue
	start = 0
	end = len(cut_positions[node])-1
	if start < end and cut_positions[node][start] < min_dist_from_end: start += 1
	if start < end and cut_positions[node][end] > nodelens[node]-1-min_dist_from_end: end -= 1
	if "<" + node in max_overlap:
		while start < end and cut_positions[node][start] <= max_overlap["<" + node]:
			sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][start]) + " due to overlap" + "\n")
			start += 1
	if ">" + node in max_overlap:
		while start < end and cut_positions[node][end] >= nodelens[node] - max_overlap[">" + node]:
			sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][end]) + " due to overlap" + "\n")
			end -= 1
	if end == start:
		cut_positions[node] = []
	else:
		cut_positions[node] = cut_positions[node][start:end+1]
		sys.stderr.write("cut sites for node " + node + " : " + " ".join(str(s) for s in cut_positions[node]) + "\n")

with open(out_mapping_file, "w") as f:
	for node in cut_positions:
		if len(cut_positions[node]) > 0:
			for i in range(0, len(cut_positions[node])):
				last_cut = 0
				if i > 0: last_cut = cut_positions[node][i-1]
				f.write(node + "_cut" + str(i) + "\t" + ">" + node + ":" + str(last_cut) + ":" + str(nodelens[node] - cut_positions[node][i]) + "\n")
			f.write(node + "_cut" + str(len(cut_positions[node])) + "\t" + ">" + node + ":" + str(cut_positions[node][-1]) + ":" + "0" + "\n")

with open(in_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			if parts[1] not in cut_positions or len(cut_positions[parts[1]]) == 0:
				print(l.strip())
			else:
				for i in range(0, len(cut_positions[parts[1]])):
					last_cut = 0
					if i > 0: last_cut = cut_positions[parts[1]][i-1]
					node_name = parts[1] + "_cut" + str(i)
					node_seq = parts[2][last_cut:-(nodelens[parts[1]] - cut_positions[parts[1]][i])]
					print("S\t" + node_name + "\t" + node_seq + "\t" + "\t".join(parts[3:]))
					print("L\t" + node_name + "\t+\t" + parts[1] + "_cut" + str(i+1) + "\t+\t0M")
				node_name = parts[1] + "_cut" + str(len(cut_positions[parts[1]]))
				node_seq = parts[2][cut_positions[parts[1]][-1]:]
				print("S\t" + node_name + "\t" + node_seq + "\t" + "\t".join(parts[3:]))
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode[1:] in cut_positions and len(cut_positions[fromnode[1:]]) > 0:
				if fromnode[0] == "<":
					fromnode = fromnode + "_cut0"
				else:
					fromnode = fromnode + "_cut" + str(len(cut_positions[fromnode[1:]]))
			if tonode[1:] in cut_positions and len(cut_positions[tonode[1:]]) > 0:
				if tonode[0] == ">":
					tonode = tonode + "_cut0"
				else:
					tonode = tonode + "_cut" + str(len(cut_positions[tonode[1:]]))
			parts[1] = fromnode[1:]
			parts[3] = tonode[1:]
			print("\t".join(parts))

with open(in_alns_file) as f:
	with open(out_cut_alns_file, "w") as f2:
		for l in f:
			parts = l.strip().split("\t")
			path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split("\t")
			realpath = []
			for node in path:
				if node[1:] not in cut_positions or len(cut_positions[node[1:]]) == 0:
					realpath.append(node)
				else:
					if node[0] == ">":
						for i in range(0, len(cut_positions[node[1:]])+1):
							realpath.append(node + "_cut" + str(i))
					else:
						for i in range(len(cut_positions[node[1:]]), -1, -1):
							realpath.append(node + "_cut" + str(i))
			parts[5] = "".join(realpath)
			f2.write("\t".join(parts) + "\n")
