#!/usr/bin/env python

import sys
import graph_functions as gf

in_graph_file = sys.argv[1]
in_paths_file = sys.argv[2]
gapname = sys.argv[3]
out_mapping_file = sys.argv[4]
# gfa to stdout

def parse_current_alns(tip_support, alns):
	if len(alns) < 2: return
	alns.sort(key=lambda x: x[0])
	for i in range(1, len(alns)):
		if alns[i-1][5] == 0 and (alns[i-1][3], alns[i][2]) in valid_tips:
			key = (alns[i-1][3], alns[i][2])
			if key not in tip_support: tip_support[key] = []
			tip_support[key].append((alns[i][4], alns[i][0] - alns[i-1][1]))
		elif alns[i][4] == 0 and (gf.revnode(alns[i][2]), gf.revnode(alns[i-1][3])) in valid_tips:
			key = (gf.revnode(alns[i][2]), gf.revnode(alns[i-1][3]))
			if key not in tip_support: tip_support[key] = []
			tip_support[key].append((alns[i-1][5], alns[i][0] - alns[i-1][1]))


node_seqs = {}
edges = {}
edge_overlaps = {}
max_overlap = {}
with open(in_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			node_seqs[parts[1]] = parts[2]
		elif parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
			edges[gf.revnode(tonode)].add(gf.revnode(fromnode))
			overlap = int(parts[5][:-1])
			edge_overlaps[gf.canon(fromnode, tonode)] = overlap
			if fromnode not in max_overlap: max_overlap[fromnode] = overlap
			if gf.revnode(tonode) not in max_overlap: max_overlap[gf.revnode(tonode)] = overlap
			max_overlap[fromnode] = max(overlap, max_overlap[fromnode])
			max_overlap[gf.revnode(tonode)] = max(overlap, max_overlap[gf.revnode(tonode)])

valid_tips = set()

for n in edges:
	if gf.revnode(n) in edges: continue
	if len(edges[n]) != 1: continue
	prev_hom = gf.revnode(gf.getone(edges[n]))
	if gf.revnode(prev_hom) not in edges or len(edges[gf.revnode(prev_hom)]) != 2: continue
	assert prev_hom in edges
	if len(edges[prev_hom]) != 2: continue
	other = None
	for edge in edges[prev_hom]:
		if edge == gf.revnode(n): continue
		other = edge
	if other is None: continue
	assert gf.revnode(other) in edges
	if len(edges[gf.revnode(other)]) != 1: continue
	valid_tips.add((gf.revnode(n), other))

node_count_in_tips = {}
for tip in valid_tips:
	if tip[1][1:] not in node_count_in_tips: node_count_in_tips[tip[1][1:]] = 0
	node_count_in_tips[tip[1][1:]] += 1
filtered_valid_tips = set()
for tip in valid_tips:
	if node_count_in_tips[tip[1][1:]] == 1: filtered_valid_tips.add(tip)
valid_tips = filtered_valid_tips

tip_support = {}

current_alns = []
current_read = ""
with open(in_paths_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] != current_read:
			parse_current_alns(tip_support, current_alns)
			current_alns = []
			current_read = parts[0]
		path = parts[5].replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		current_alns.append((int(parts[2]), int(parts[3]), path[0], path[-1], int(parts[7]), int(parts[6]) - int(parts[8])))
parse_current_alns(tip_support, current_alns)

node_cuts = {}
node_trims = {}
num_fixed = 0
extra_edges = []
extra_nocut_edges = []
for key in tip_support:
	if len(tip_support[key]) < 1: continue
	if (key[1], key[0]) in tip_support and len(tip_support[(key[1], key[0])]) > len(tip_support[key]):
		sys.stderr.write("conflicting fixes at " + key[0] + " " + key[1] + ", skipping lower support" + "\n")
		continue
	if (key[1], key[0]) in tip_support and len(tip_support[(key[1], key[0])]) == len(tip_support[key]):
		sys.stderr.write("ambiguous fix at " + key[0] + " " + key[1] + ", skipping" + "\n")
		continue
	wanted_cut_pos = -1
	wanted_gap_length = 0
	for t in tip_support[key]:
		if t[0] > wanted_cut_pos:
			wanted_cut_pos = t[0]
			wanted_gap_length = t[1]
	assert wanted_cut_pos >= 0
	if abs(wanted_gap_length) >= len(node_seqs[key[0][1:]]) or abs(wanted_gap_length) >= len(node_seqs[key[1][1:]]):
		sys.stderr.write("can't fix " + key[0] + " " + key[1] + " due to overlap containing node (wanted " + str(wanted_gap_length) + ", node lengths " + str(len(node_seqs[key[0][1:]])) + ", " + str(len(node_seqs[key[1][1:]])) + ")")
		continue
	if wanted_cut_pos == 0:
		sys.stderr.write("mend " + key[0] + " " + key[1] + "\n")
		extra_nocut_edges.append((key[0], key[1], wanted_gap_length))
		num_fixed += 1
		continue
	max_ov = 0
	if gf.revnode(key[1]) in max_overlap: max_ov = max_overlap[gf.revnode(key[1])]
	max_rev_ov = 0
	if key[1] in max_overlap: max_rev_ov = max_overlap[key[1]]
	if wanted_cut_pos <= max_ov:
		sys.stderr.write("can't fix " + key[0] + " " + key[1] + " due to overlap (wanted " + str(wanted_cut_pos) + ", overlap " + str(max_ov) + ")" + "\n")
		continue
	if wanted_cut_pos >= len(node_seqs[key[1][1:]]) - max_rev_ov:
		sys.stderr.write("can't fix " + key[0] + " " + key[1] + " due to overlap (wanted " + str(len(node_seqs[key[1][1:]]) - wanted_cut_pos) + ", overlap " + str(max_rev_ov) + ")" + "\n")
		continue
	if wanted_gap_length < 0:
		if -wanted_gap_length >= len(node_seqs[key[1][1:]]) - wanted_cut_pos:
			wanted_trim = -wanted_gap_length - (len(node_seqs[key[1][1:]]) - wanted_cut_pos) + 1
			assert wanted_trim >= 1
			sys.stderr.write("trim back " + str(key[0]) + " by " + str(wanted_trim) + "\n")
			assert gf.revnode(key[0]) in max_overlap
			assert wanted_trim < len(node_seqs[key[0][1:]])
			if len(node_seqs[key[0][1:]]) - wanted_trim <= max_overlap[gf.revnode(key[0])]:
				sys.stderr.write("can't trim due to overlap, skipping")
				continue
			node_trims[key[0]] = wanted_trim
			wanted_gap_length += wanted_trim
	sys.stderr.write("mend " + key[0] + " " + key[1] + "\n")
	assert key[1] not in node_cuts
	assert gf.revnode(key[1]) not in node_cuts
	node_cuts[key[1]] = wanted_cut_pos
	extra_edges.append((key[0], key[1], wanted_gap_length))
	num_fixed += 1

for node in node_seqs:
	if ">" + node not in node_cuts and "<" + node not in node_cuts and ">" + node not in node_trims and "<" + node not in node_trims:
		print("S\t" + node + "\t" + node_seqs[node])
	if (">" + node) in node_trims:
		assert ">" + node not in node_cuts
		assert "<" + node not in node_cuts
		assert "<" + node not in node_trims
		print("S\t" + node + "_trim" + "\t" + node_seqs[node][:-node_trims[">" + node]])
	if ("<" + node) in node_trims:
		assert ">" + node not in node_cuts
		assert "<" + node not in node_cuts
		assert ">" + node not in node_trims
		print("S\t" + node + "_trim" + "\t" + node_seqs[node][node_trims["<" + node]:])
	if (">" + node) in node_cuts:
		assert ">" + node not in node_trims
		assert "<" + node not in node_trims
		assert "<" + node not in node_cuts
		cut_pos = node_cuts[">" + node]
		print("S\t" + node + "_beg" + "\t" + node_seqs[node][:cut_pos])
		print("L\t" + node + "_beg" + "\t" + "+" + "\t" + node + "_end" + "\t" + "+" + "\t" + "0M")
		print("S\t" + node + "_end" + "\t" + node_seqs[node][cut_pos:])
	if ("<" + node) in node_cuts:
		assert ">" + node not in node_trims
		assert "<" + node not in node_trims
		assert ">" + node not in node_cuts
		cut_pos = len(node_seqs[node]) - node_cuts["<" + node]
		print("S\t" + node + "_beg" + "\t" + node_seqs[node][:cut_pos])
		print("L\t" + node + "_beg" + "\t" + "+" + "\t" + node + "_end" + "\t" + "+" + "\t" + "0M")
		print("S\t" + node + "_end" + "\t" + node_seqs[node][cut_pos:])

for edge in edges:
	for edge2 in edges[edge]:
		fromnode = edge
		tonode = edge2
		overlap=edge_overlaps[gf.canon(edge, edge2)]
		if ">" + fromnode[1:] in node_trims or "<" + fromnode[1:] in node_trims:
			if ">" + fromnode[1:] in node_trims and ">" + fromnode[1:] in gf.canon(edge, edge2): overlap -= node_trims[">" + fromnode[1:]]
			if "<" + fromnode[1:] in node_trims and "<" + fromnode[1:] in gf.canon(edge, edge2): overlap -= node_trims["<" + fromnode[1:]]

			fromnode = fromnode + "_trim"
		if ">" + tonode[1:] in node_trims or "<" + tonode[1:] in node_trims:
			if ">" + tonode[1:] in node_trims and ">" + tonode[1:] in gf.canon(edge, edge2): overlap -= node_trims[">" + tonode[1:]]
			if "<" + tonode[1:] in node_trims and "<" + tonode[1:] in gf.canon(edge, edge2): overlap -= node_trims["<" + tonode[1:]]

			tonode = tonode + "_trim"
		if fromnode[0] == ">" and (edge in node_cuts or gf.revnode(edge) in node_cuts):
			fromnode = fromnode + "_end"
		elif fromnode[0] == "<" and (edge in node_cuts or gf.revnode(edge) in node_cuts):
			fromnode = fromnode + "_beg"
		if tonode[0] == ">" and (edge2 in node_cuts or gf.revnode(edge2) in node_cuts):
			tonode = tonode + "_beg"
		elif tonode[0] == "<" and (edge2 in node_cuts or gf.revnode(edge2) in node_cuts):
			tonode = tonode + "_end"
		print("L\t" + fromnode[1:] + "\t" + ("+" if fromnode[0] == ">" else "-") + "\t" + tonode[1:] + "\t" + ("+" if tonode[0] == ">" else "-") + "\t" + str(overlap) + "M")

next_fake_node_id = 0
for edge in extra_edges:
	assert gf.revnode(edge[0]) not in node_trims
	fromprint = ""
	if edge[0] in node_cuts or gf.revnode(edge[0]) in node_cuts:
		fromprint = "_end" if edge[0][0] == ">" else "_beg"
		assert edge[0] not in node_trims
	assert edge[1] in node_cuts
	assert gf.revnode(edge[1]) not in node_cuts
	gap_len = edge[2]
	trimprint = ""
	if edge[0] in node_trims: trimprint = "_trim"
	if gap_len <= 0:
		print("L\t" + edge[0][1:] + fromprint + trimprint + "\t" + ("+" if edge[0][0] == ">" else "-") + "\t" + edge[1][1:] + ("_end" if edge[1][0] == ">" else "_beg") + "\t" + ("+" if edge[1][0] == ">" else "-") + "\t" + str(-gap_len) + "M")
	else:
		assert trimprint == ""
		assert edge[0] not in node_trims
		fake_node = gapname + "_" + str(next_fake_node_id)
		next_fake_node_id += 1
		print("L\t" + edge[0][1:] + fromprint + "\t" + ("+" if edge[0][0] == ">" else "-") + "\t" + fake_node + "\t" + "+" + "\t" + "0M")
		print("S\t" + fake_node + "\t" + "N"*gap_len)
		print("L\t" + fake_node + "\t" + "+" + "\t" + edge[1][1:] + ("_end" if edge[1][0] == ">" else "_beg") + "\t" + ("+" if edge[1][0] == ">" else "-") + "\t" + "0M")

for edge in extra_nocut_edges:
	assert edge[0] not in node_cuts
	assert gf.revnode(edge[0]) not in node_cuts
	assert edge[1] not in node_cuts
	assert gf.revnode(edge[1]) not in node_cuts
	gap_len = edge[2]
	trimprint = ""
	if edge[0] in node_trims: trimprint = "_trim"
	if gap_len <= 0:
		print("L\t" + edge[0][1:] + trimprint + "\t" + ("+" if edge[0][0] == ">" else "-") + "\t" + edge[1][1:] + "\t" + ("+" if edge[1][0] == ">" else "-") + "\t" + str(-gap_len) + "M")
	else:
		fake_node = gapname + "_" + str(next_fake_node_id)
		next_fake_node_id += 1
		print("L\t" + edge[0][1:] + "\t" + ("+" if edge[0][0] == ">" else "-") + "\t" + fake_node + "\t" + "+" + "\t" + "0M")
		print("S\t" + fake_node + "\t" + "N"*gap_len)
		print("L\t" + fake_node + "\t" + "+" + "\t" + edge[1][1:] + "\t" + ("+" if edge[1][0] == ">" else "-") + "\t" + "0M")

with open(out_mapping_file, "w") as f:
	for cut in node_cuts:
		node = cut[1:]
		if cut[0] == ">":
			cut_pos = node_cuts[">" + node]
		else:
			assert cut[0] == "<"
			cut_pos = len(node_seqs[node]) - node_cuts["<" + node]
		f.write(node + "_beg" + "\t" + ">" + node + ":0:" + str(len(node_seqs[node]) - cut_pos) + "\n")
		f.write(node + "_end" + "\t" + ">" + node + ":" + str(cut_pos) + ":0" + "\n")
	for trim in node_trims:
		node = trim[1:]
		if trim[0] == ">":
			f.write(node + "_trim" + "\t" + ">" + node + ":0:" + str(node_trims[trim]) + "\n")
		else:
			f.write(node + "_trim" + "\t" + ">" + node + ":" + str(node_trims[trim]) + ":0" + "\n")

