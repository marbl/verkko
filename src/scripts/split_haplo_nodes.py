#!/usr/bin/env python

import sys

in_gfa = sys.argv[1]
in_matches = sys.argv[2]
in_node_alns = sys.argv[3]
in_read_alns = sys.argv[4]
in_distance_cuts = sys.argv[5]
out_split_gfa = sys.argv[6]
out_relevant_read_names = sys.argv[7]


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
				cut_pos = int(parts[6]) - int(parts[8])
		else:
			if parts[4] == "+":
				cut_pos = int(parts[7])
			else:
				cut_pos = int(parts[6]) - int(parts[7])
		if cut_direction[parts[0]]:
			if parts[0] not in best_cut or int(parts[1]) - int(parts[3]) < best_cut[parts[0]][0]:
				best_cut[parts[0]] = (int(parts[1]) - int(parts[3]), parts[5], cut_pos, parts[4])
		else:
			if parts[0] not in best_cut or int(parts[2]) < best_cut[parts[0]][0]:
				best_cut[parts[0]] = (int(parts[2]), parts[5], cut_pos, parts[4])

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

cuts = []
for tip_node in best_cut:
	cut_offset = best_cut[tip_node][0]
	cut_node = best_cut[tip_node][1]
	cut_pos = best_cut[tip_node][2]
	node_direction = best_cut[tip_node][3]
	cuts.append((tip_node, cut_offset, cut_node, cut_pos, node_direction))

with open(in_distance_cuts) as f:
	for l in f:
		parts = l.strip().split('\t')
		tip_node = parts[0][1:]
		assert (cut_direction[tip_node]) == (parts[0][0] == ">")
		cut_offset = 0
		cut_node = parts[1][1:]
		cut_pos = int(parts[2])
		if parts[1][0] == "<": cut_pos = nodelens[parts[1][1:]] - cut_pos
		node_direction = ("+" if (parts[1][0] == parts[0][0]) else "-")
		cuts.append((tip_node, cut_offset, cut_node, cut_pos, node_direction))

cut_poses = {}
for cut_triplet in cuts:
	tip_node = cut_triplet[0]
	cut_offset = cut_triplet[1]
	cut_node = cut_triplet[2]
	cut_pos = cut_triplet[3]
	node_direction = cut_triplet[4]
	tip_direction = cut_direction[tip_node]
	cut_node_direction = node_direction == "+"
	if tip_direction == cut_node_direction:
		cut_pos += cut_offset
	else:
		cut_pos -= cut_offset
	assert ">" + cut_node in max_edge_overlaps
	assert "<" + cut_node in max_edge_overlaps
	try_backtrace = 0
	if max_edge_overlaps["<" + cut_node] >= cut_pos:
		try_backtrace = max_edge_overlaps["<" + cut_node] - cut_pos + 1
		assert try_backtrace > 0
		sys.stderr.write("Shifted cut site " + cut_node + " " + str(cut_pos) + " to " + str(cut_pos + try_backtrace) + " due to overlaps" + "\n")
		cut_pos = cut_pos + try_backtrace
	if nodelens[cut_node] -  max_edge_overlaps[">" + cut_node] <= cut_pos:
		try_backtrace = -(cut_pos - (nodelens[cut_node] - max_edge_overlaps[">" + cut_node]) + 1)
		assert try_backtrace < 0
		sys.stderr.write("Shifted cut site " + cut_node + " " + str(cut_pos) + " to " + str(cut_pos + try_backtrace) + " due to overlaps" + "\n")
		cut_pos = cut_pos + try_backtrace
	if max_edge_overlaps["<" + cut_node] >= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	if nodelens[cut_node] - max_edge_overlaps[">" + cut_node] <= cut_pos:
		sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to edge overlap" + "\n")
		continue
	if try_backtrace != 0:
		if tip_direction:
			tip_cut_position = nodelens[tip_node]
		else:
			tip_cut_position = 0
		if cut_node_direction:
			tip_cut_position += try_backtrace
		else:
			tip_cut_position -= try_backtrace
		if max_edge_overlaps["<" + tip_node] >= tip_cut_position:
			sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to tip edge overlap (" + str(tip_cut_position) + ")" + "\n")
			continue
		if nodelens[tip_node] - max_edge_overlaps[">" + tip_node] <= tip_cut_position:
			sys.stderr.write("Discarded cut site " + cut_node + " " + str(cut_pos) + " due to tip edge overlap (" + str(tip_cut_position) + ")" + "\n")
			continue
		if tip_node not in cut_poses: cut_poses[tip_node] = []
		cut_poses[tip_node].append(tip_cut_position)
		sys.stderr.write("Cut " + tip_node + " " + str(tip_cut_position) + "\n")
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
				f2.write("S\t" + parts[1] + "_hapcutfirst\t" + parts[2][0:cuts[0]] + "\n")
				for i in range(1, len(cuts)):
					f2.write("S\t" + parts[1] + "_hapcut" + str(i) + "\t" + parts[2][cuts[i-1]:cuts[i]] + "\n")
				f2.write("S\t" + parts[1] + "_hapcutlast\t" + parts[2][cuts[-1]:] + "\n")
				continue
			if parts[0] == "L" and (parts[1] in cut_poses or parts[3] in cut_poses):
				if parts[1] in cut_poses:
					if parts[2] == "+":
						parts[1] = parts[1] + "_hapcutlast"
					else:
						parts[1] = parts[1] + "_hapcutfirst"
				if parts[3] in cut_poses:
					if parts[4] == "+":
						parts[3] = parts[3] + "_hapcutfirst"
					else:
						parts[3] = parts[3] + "_hapcutlast"
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
