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
		if readname not in alns_per_read: alns_per_read[readname] = 0
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
		start_node = path[0]
		if path[0][0] == "<": start_node_offset = nodelens[start_node[1:]] - start_node_offset - 1
		end_node = revnode(path[-1])
		end_node_offset = int(parts[6]) - int(parts[8])
		if path[-1][0] == ">": end_node_offset = nodelens[end_node[1:]] - end_node_offset - 1
		if readname not in read_aln_positions: read_aln_positions[readname] = []
		read_aln_positions[readname].append((alnstart, alnend, start_node, start_node_offset, end_node, end_node_offset))

aln_connect_positions = {}
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
		if prev_node > this_node or (prev_node == this_node and prev_node_offset > this_node_offset):
			(prev_node, this_node) = (this_node, prev_node)
			(prev_node_offset, this_node_offset) = (this_node_offset, prev_node_offset)
		key = (prev_node, this_node)
		if key not in aln_connect_positions: aln_connect_positions[key] = []
		aln_connect_positions[key].append((prev_node_offset, this_node_offset))

cut_positions = {}
for key in aln_connect_positions:
	poses = aln_connect_positions[key]
	poses.sort(key=lambda x: x[1])
	poses.sort(key=lambda x: x[0])
	cluster_origin = []
	cluster_poses = []
	for i in range(0, len(poses)):
		cluster_origin.append(i)
		cluster_poses.append([poses[i]])
	for i in range(0, len(poses)):
		for j in range(i-1, -1, -1):
			if abs(poses[i][0] - poses[j][0]) <= cut_length_threshold and abs(poses[i][1] - poses[j][1]) <= cut_length_threshold:
				cluster_origin[i] = j
				break
	for i in range(len(poses)-1, -1, -1):
		if cluster_origin[i] == i:
			if len(cluster_poses[i]) >= cut_coverage_threshold:
				start_poses = [x[0] for x in cluster_poses[i]]
				end_poses = [x[1] for x in cluster_poses[i]]
				if key[0][1:] not in cut_positions: cut_positions[key[0][1:]] = []
				if key[1][1:] not in cut_positions: cut_positions[key[1][1:]] = []
				cut_positions[key[0][1:]].append(start_poses[len(start_poses)//2])
				cut_positions[key[1][1:]].append(end_poses[len(end_poses)//2])
		else:
			assert cluster_origin[i] < i
			cluster_poses[cluster_origin[i]] += cluster_poses[i]

for node in cut_positions:
	if len(cut_positions[node]) < 1: continue
	cut_positions[node].sort()
	start = 0
	end = len(cut_positions[node])-1
	while start <= end and cut_positions[node][start] < min_dist_from_end:
		sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][start]) + " due to too close to start" + "\n")
		start += 1
	while start <= end and cut_positions[node][end] > nodelens[node]-1-min_dist_from_end:
		sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][end]) + " due to too close to end" + "\n")
		end -= 1
	if "<" + node in max_overlap:
		while start <= end and cut_positions[node][start] <= max_overlap["<" + node]:
			sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][start]) + " due to overlap" + "\n")
			start += 1
	if ">" + node in max_overlap:
		while start <= end and cut_positions[node][end] >= nodelens[node] - max_overlap[">" + node]:
			sys.stderr.write("discarded cut site at node " + str(node) + " at position " + str(cut_positions[node][end]) + " due to overlap" + "\n")
			end -= 1
	if end < start:
		cut_positions[node] = []
	else:
		chosen_positions = []
		for i in range(start, end+1):
			if len(chosen_positions) >= 1 and cut_positions[node][i] == chosen_positions[-1]:
				continue
			assert len(chosen_positions) == 0 or cut_positions[node][i] > chosen_positions[-1]
			chosen_positions.append(cut_positions[node][i])
		cut_positions[node] = chosen_positions
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
