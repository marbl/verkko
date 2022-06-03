#!/usr/bin/env python

import sys

in_gfa = sys.argv[1]
in_read_alns = sys.argv[2]
in_distance_cuts = sys.argv[3]
out_split_gfa = sys.argv[4]
out_relevant_read_names = sys.argv[5]
out_split_mapping = sys.argv[6]


def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

nodelens = {}
max_edge_overlaps = {}
with open(in_gfa) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			nodelens[parts[1]] = len(parts[2])
			if ">" + parts[1] not in max_edge_overlaps: max_edge_overlaps[">" + parts[1]] = 0
			if "<" + parts[1] not in max_edge_overlaps: max_edge_overlaps["<" + parts[1]] = 0
		if parts[0] == "L":
			overlap = int(parts[5][:-1])
			fromtip = (">" if parts[2] == "+" else "<") + parts[1]
			totip = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromtip not in max_edge_overlaps: max_edge_overlaps[fromtip] = 0
			if totip not in max_edge_overlaps: max_edge_overlaps[totip] = 0
			max_edge_overlaps[fromtip] = max(max_edge_overlaps[fromtip], overlap)
			max_edge_overlaps[totip] = max(max_edge_overlaps[totip], overlap)

relevant_nodes = set()
cuts = []
with open(in_distance_cuts) as f:
	for l in f:
		parts = l.strip().split('\t')
		tip_node = parts[0]
		tip_cut_pos = nodelens[parts[0][1:]]
		cut_node = parts[1]
		cut_pos = int(parts[2])
		cuts.append((tip_node, tip_cut_pos, cut_node, cut_pos))
		relevant_nodes.add(tip_node[1:])
		relevant_nodes.add(cut_node[1:])

cut_poses = {}
for cut_triplet in cuts:
	tip_node = cut_triplet[0]
	tip_cut_pos = cut_triplet[1]
	cut_node = cut_triplet[2]
	cut_pos = cut_triplet[3]
	try_backtrace = 0
	if max_edge_overlaps[revnode(cut_node)] >= cut_pos:
		try_backtrace = max_edge_overlaps[revnode(cut_node)] - cut_pos + 1
		assert try_backtrace > 0
		sys.stderr.write("Shifted cut site " + cut_node + " " + str(cut_pos) + " to " + str(cut_pos + try_backtrace) + " due to overlaps" + "\n")
	if nodelens[cut_node[1:]] - max_edge_overlaps[cut_node] <= cut_pos:
		try_backtrace = -(cut_pos - (nodelens[cut_node[1:]] - max_edge_overlaps[cut_node]) + 1)
		assert try_backtrace < 0
		sys.stderr.write("Shifted cut site " + cut_node + " " + str(cut_pos) + " to " + str(cut_pos + try_backtrace) + " due to overlaps" + "\n")
	cut_pos = cut_pos + try_backtrace
	tip_cut_pos += try_backtrace
	if max_edge_overlaps[revnode(cut_node)] >= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	if nodelens[cut_node[1:]] - max_edge_overlaps[cut_node] <= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	assert max_edge_overlaps[tip_node] == 0
	if tip_cut_pos != nodelens[tip_node[1:]]:
		if tip_cut_pos > nodelens[tip_node[1:]]:
			sys.stderr.write("Discarded cut site " + cut_node + " " + str(tip_cut_pos) + " due to outside tip" + "\n")
			continue
		if tip_cut_pos <= max_edge_overlaps[revnode(tip_node)]:
			sys.stderr.write("Discarded cut site " + cut_node + " " + str(tip_cut_pos) + " due to outside tip" + "\n")
			continue
	assert tip_cut_pos > 0
	assert tip_cut_pos <= nodelens[tip_node[1:]]
	if tip_cut_pos != nodelens[tip_node[1:]]:
		if tip_node[1:] not in cut_poses: cut_poses[tip_node[1:]] = []
		tip_cut_position = tip_cut_pos
		if tip_node[0] == "<": tip_cut_position = nodelens[tip_node[1:]] - tip_cut_position
		cut_poses[tip_node[1:]].append(tip_cut_position)
		sys.stderr.write("Cut " + tip_node[1:] + " " + str(tip_cut_position) + "\n")
	assert cut_pos > 0
	assert cut_pos < nodelens[cut_node[1:]]
	if cut_node[1:] not in cut_poses: cut_poses[cut_node[1:]] = []
	if cut_node[0] == "<": cut_pos = nodelens[cut_node[1:]] - cut_pos
	cut_poses[cut_node[1:]].append(cut_pos)
	sys.stderr.write("Cut " + cut_node[1:] + " " + str(cut_pos) + "\n")

mappings = []
cuts_per_node = {}
for node in cut_poses:
	cuts_per_node[node] = len(set(cut_poses[node]))
with open(out_split_gfa, "w") as f2:
	with open(in_gfa) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == "S" and parts[1] in cut_poses:
				nodelen = len(parts[2])
				cuts = list(set(cut_poses[parts[1]]))
				cuts.sort()
				assert parts[1] in cuts_per_node
				assert len(cuts) == cuts_per_node[parts[1]]
				assert len(cuts) >= 1
				f2.write("S\t" + parts[1] + "_hapcut0\t" + parts[2][0:cuts[0]] + "\n")
				mappings.append((parts[1] + "_hapcut0", parts[1], 0, nodelens[parts[1]] - cuts[0]))
				for i in range(1, len(cuts)):
					f2.write("S\t" + parts[1] + "_hapcut" + str(i) + "\t" + parts[2][cuts[i-1]:cuts[i]] + "\n")
					mappings.append((parts[1] + "_hapcut" + str(i), parts[1], cuts[i-1], nodelens[parts[1]] - cuts[i]))
				f2.write("S\t" + parts[1] + "_hapcut" + str(len(cuts)) + "\t" + parts[2][cuts[-1]:] + "\n")
				mappings.append((parts[1] + "_hapcut" + str(len(cuts)), parts[1], cuts[-1], 0))
				continue
			if parts[0] == "L" and (parts[1] in cut_poses or parts[3] in cut_poses):
				if parts[1] in cut_poses:
					if parts[2] == "+":
						parts[1] = parts[1] + "_hapcut" + str(cuts_per_node[parts[1]])
					else:
						parts[1] = parts[1] + "_hapcut0"
				if parts[3] in cut_poses:
					if parts[4] == "+":
						parts[3] = parts[3] + "_hapcut0"
					else:
						parts[3] = parts[3] + "_hapcut" + str(cuts_per_node[parts[3]])
				f2.write("\t".join(parts) + "\n")
				continue
			else:
				f2.write(l.strip() + "\n")

with open(out_relevant_read_names, "w") as f2:
	with open(in_read_alns) as f:
		for l in f:
			parts = l.strip().split('\t')
			path = parts[5].replace('<', '\t').replace('>', '\t').strip().split('\t')
			for node in path:
				if node in relevant_nodes:
					f2.write(parts[0].split(' ')[0] + "\n")

with open(out_split_mapping, "w") as f:
	for t in mappings:
		f.write(t[0] + "\t>" + t[1] + ":" + str(t[2]) + ":" + str(t[3]) + "\n")
