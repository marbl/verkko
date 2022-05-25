#!/usr/bin/env python

import sys

in_gfa = sys.argv[1]
in_matches = sys.argv[2]
in_node_alns = sys.argv[3]
in_read_alns = sys.argv[4]
out_split_gfa = sys.argv[5]
out_relevant_read_names = sys.argv[6]


relevant_nodes = set()
match_pairs = set()
cut_direction = {}
with open(in_matches) as f:
	for l in f:
		parts = l.strip().split('\t')
		relevant_nodes.add(parts[0][1:])
		for i in parts[1:]:
			relevant_nodes.add(i[1:])
			match_pairs.add((parts[0][1:], i[1:], parts[0][0] == i[0]))
		assert parts[0][1:] not in cut_direction
		cut_direction[parts[0][1:]] = (parts[0][0] == ">")

best_cut = {}
with open(in_node_alns) as f:
	for l in f:
		parts = l.strip().split('\t')
		if (parts[0], parts[5], parts[4] == "+") not in match_pairs: continue
		if cut_direction[parts[0]]:
			if parts[4] == "+":
				cut_pos = int(parts[8])
			else:
				cut_pos = int(parts[6]) - int(parts[7])
		else:
			if parts[4] == "-":
				cut_pos = int(parts[6]) - int(parts[8])
			else:
				cut_pos = int(parts[7])
		if cut_direction[parts[0]]:
			if parts[0] not in best_cut or int(parts[3]) > best_cut[parts[0]][0]:
				best_cut[parts[0]] = (int(parts[3]), parts[5], cut_pos)
		else:
			if parts[0] not in best_cut or int(parts[2]) < best_cut[parts[0]][0]:
				best_cut[parts[0]] = (int(parts[2]), parts[5], cut_pos)

nodelens = {}
max_edge_overlaps = {}
with open(in_gfa) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			nodelens[parts[1]] = len(parts[2])
		if parts[0] == "L":
			overlap = int(parts[5][:-1])
			fromtip = (">" if parts[2] == "+" else "<") + parts[1]
			totip = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromtip not in max_edge_overlaps: max_edge_overlaps[fromtip] = 0
			if totip not in max_edge_overlaps: max_edge_overlaps[totip] = 0
			max_edge_overlaps[fromtip] = max(max_edge_overlaps[fromtip], overlap)
			max_edge_overlaps[totip] = max(max_edge_overlaps[totip], overlap)

cut_poses = {}
for tip_node in best_cut:
	cut_node = best_cut[tip_node][1]
	cut_pos = best_cut[tip_node][2]
	if "<" + cut_node in max_edge_overlaps and max_edge_overlaps["<" + cut_node] >= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	if ">" + cut_node in max_edge_overlaps and nodelens[cut_node] - max_edge_overlaps[">" + cut_node] <= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	if cut_node not in cut_poses: cut_poses[cut_node] = []
	cut_poses[cut_node].append(cut_pos)
	sys.stderr.write("Cut " + cut_node + " " + str(cut_pos) + "\n")

with open(out_split_gfa, "w") as f2:
	with open(in_gfa) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == "S" and parts[1] in cut_poses:
				nodelen = len(parts[2])
				cuts = cut_poses[parts[1]]
				cuts.sort()
				assert len(cuts) >= 1
				f2.write("S\t" + parts[1] + "_first\t" + parts[2][0:cuts[0]] + "\n")
				for i in range(1, len(cuts)-1):
					f2.write("S\t" + parts[1] + "_cut" + str(i) + "\t" + parts[2][cuts[i-1]:cuts[i]] + "\n")
				f2.write("S\t" + parts[1] + "_last\t" + parts[2][cuts[-1]:] + "\n")
				continue
			if parts[0] == "L" and (parts[1] in cut_poses or parts[3] in cut_poses):
				if parts[1] in cut_poses:
					if parts[2] == "+":
						parts[1] = parts[1] + "_last"
					else:
						parts[1] = parts[1] + "_first"
				if parts[3] in cut_poses:
					if parts[4] == "+":
						parts[3] = parts[3] + "_first"
					else:
						parts[3] = parts[3] + "_last"
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
