#!/usr/bin/env python

import sys
import networkx as nx
import graph_functions as gf

graph_file = sys.argv[1]
input_aln_file = sys.argv[2]
min_gap_coverage = int(sys.argv[3])
max_end_clip = int(sys.argv[4])
nongap_aln_file_out = sys.argv[5]
gap_aln_file_out = sys.argv[6]
prefix = sys.argv[7]
allow_nontips = (True if sys.argv[8] == "y" else False)
allow_onetip = (True if sys.argv[8] == "o" else False)
max_read_clip = 0.05
# graph to stdout

not_tips = set()
node_seqs = {}
edge_overlaps = {}

#not_tips = no outgoing nodes, assymetric!

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			node_seqs[parts[1]] = parts[2]
		if parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
#this is not to_node but reverse(to_node)
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			edge_overlaps[gf.canontip(fromnode, tonode)] = int(parts[5][:-1])
			if parts[1] != parts[3]: # we will skip self-edges for this consideration, if you have a loop at a gap, try to patch it anyway
				not_tips.add(fromnode)
				not_tips.add(tonode)
		print(l.strip())

if allow_nontips: not_tips = set()


G = nx.Graph()
gf.load_indirect_graph(graph_file, G)
tangles = gf.nodes_in_tangles(G, 30000, 100)
or_tangles = set()
for node in tangles:
	for pref in (">", "<"):
		or_tangles.add(pref + node)

gaps = {}

count_per_read = {}
with open(input_aln_file) as f:
	for l in f:
		readname = l.strip().split('\t')[0].split(' ')[0]
		if readname not in count_per_read: count_per_read[readname] = 0
		count_per_read[readname] += 1

#Here we save rc to start of alignment. This allow to work 
alns_per_read = {}
with open(input_aln_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0].split(' ')[0]
		if count_per_read[readname] < 2: continue
		readlen = int(parts[1])
		readstart = int(parts[2])
		readend = int(parts[3])
		path = tuple(parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t'))
		alnstart = gf.revnode(path[0])
		alnend = path[-1]
		leftclip = int(parts[7])
		rightclip = int(parts[6]) - int(parts[8])
		assert rightclip >= 0
		# skip an alignment if it's in the middle of a read and not spanning a contig
		if (readstart > int(max_read_clip*readlen) and readlen - readend > int(max_read_clip*readlen) and (leftclip > max_end_clip or rightclip > max_end_clip)): continue
		# skip an alignment if the ends are too far from the end of a contig
		if leftclip > max_end_clip and rightclip > max_end_clip: continue
		if readname not in alns_per_read: alns_per_read[readname] = []
		alns_per_read[readname].append((readstart, readend, alnstart, alnend, leftclip, rightclip, readlen, path))

used_names = list(alns_per_read)
used_names.sort()

def tangle_onetip(s_node, e_node, or_tangles, not_tips):
	if (s_node in not_tips and e_node in or_tangles) or (e_node in not_tips and s_node in or_tangles):
#		print (f"banning nodes near tangles {s_node} {e_node}")
#		print (f"{s_node in not_tips} and {e_node in or_tangles} or {e_node in not_tips} ad {s_node in or_tangles}")
		return True
	else:
		return False
		

for name in alns_per_read:
	alns_per_read[name].sort(key=lambda x: x[0])
	if len(alns_per_read[name]) < 2: continue
	alns = alns_per_read[name]
	for i in range(1, len(alns)):
		assert alns[i][0] >= alns[i-1][0]
		if alns[i][0] == alns[i-1][0]: continue
		if alns[i][1] <= alns[i-1][1]: continue
#alns[i][2] is RC to alignment, that's why it works.		
		if allow_onetip and alns[i-1][3] in not_tips and alns[i][2] in not_tips: continue	
		if allow_onetip and tangle_onetip(alns[i-1][3], alns[i][2], or_tangles, not_tips): continue 			
		if not allow_onetip and (alns[i-1][3] in not_tips or alns[i][2] in not_tips): continue
		if alns[i-1][5] > max_end_clip or alns[i][4] > max_end_clip: continue 
		gap_len = (alns[i][0] - alns[i][4]) - (alns[i-1][1] + alns[i-1][5])
		key = gf.canontip(alns[i-1][3], alns[i][2])
		# skip inserting gaps between contigs that already have edges or to ourselves
		if alns[i-1][3][1:] == alns[i][2][1:] or key in edge_overlaps: continue
		if key not in gaps: gaps[key] = []
		gaps[key].append(gap_len)
		readlen = alns[i][6]
		readstart = alns[i-1][1] + alns[i-1][5]
		readend = alns[i][0] - alns[i][4]

next_gap_id = 1
gap_lens = {}
allowed_gaps = set()
gap_names = {}
gapkeys = list(gaps)
gapkeys.sort()

for gap in gapkeys:
	if len(gaps[gap]) < min_gap_coverage: continue
	gap_len_list = list(gaps[gap])
	gap_len_list.sort()
	gap_len = gap_len_list[len(gap_len_list) // 2]
	gap_name = prefix + "-" + str(next_gap_id) + "-len-" + str(gap_len) + "-cov-" + str(len(gaps[gap]))
	gap_names[gap] = gap_name
	overlap = 0
	seq = "N"
	if gap_len < 0:
		overlap = -gap_len
		if overlap >= len(node_seqs[gap[0][1:]])-1: overlap = len(node_seqs[gap[0][1:]])-2
		if overlap >= len(node_seqs[gap[1][1:]])-1: overlap = len(node_seqs[gap[1][1:]])-2
		base_seq = node_seqs[gap[0][1:]]
		if gap[0][0] == "<": base_seq = gf.revcomp(base_seq)
		assert len(base_seq) > overlap+1
		seq = base_seq[-overlap-1:]
		base_seq = node_seqs[gap[1][1:]]
		if gap[1][0] == ">": base_seq = gf.revcomp(base_seq)
		assert len(base_seq) > overlap+1
		seq += base_seq[overlap]
	elif gap_len == 0:
		base_seq = node_seqs[gap[0][1:]]
		if gap[0][0] == "<": base_seq = gf.revcomp(base_seq)
		seq = base_seq[-1]
		base_seq = node_seqs[gap[1][1:]]
		if gap[1][0] == ">": base_seq = gf.revcomp(base_seq)
		seq += base_seq[0]
		overlap = 1
	else:
		seq = "N" * gap_len
		overlap = 0
	gap_lens[gap_name] = len(seq)
	if gap_len < 1: gap_len = 1
	next_gap_id += 1
	node_seqs[gap_name] = seq
	edge_overlaps[gf.canontip(gap[0], "<" + gap_name)] = overlap
	edge_overlaps[gf.canontip(">" + gap_name, gap[1])] = overlap
	overlap = str(overlap) + "M"
	print("S\t" + gap_name + "\t" + seq)
	print("L\t" + gap[0][1:] + "\t" + ("+" if gap[0][0] == ">" else "-") + "\t" + gap_name + "\t+\t" + overlap)
	print("L\t" + gap_name + "\t+\t" + gap[1][1:] + "\t" + ("-" if gap[1][0] == ">" else "+") + "\t" + overlap)
	allowed_gaps.add(gap)

used_gap_alns = set()

with open(gap_aln_file_out, "w") as f:
	for name in used_names:
		if len(alns_per_read[name]) < 2: continue
		alns = alns_per_read[name]
		last_used_index = None
		current_read_start = None
		current_read_end = None
		current_path = None
		current_left_clip = None
		current_right_clip = None
		for i in range(1, len(alns)):
			assert alns[i][0] >= alns[i-1][0]
			if alns[i][0] == alns[i-1][0]: continue
			if alns[i][1] <= alns[i-1][1]: continue
			if allow_onetip and alns[i-1][3] in not_tips and alns[i][2] in not_tips: continue
			if allow_onetip and tangle_onetip(alns[i-1][3], alns[i][2], or_tangles, not_tips): continue			
			if not allow_onetip and (alns[i-1][3] in not_tips or alns[i][2] in not_tips): continue
			if alns[i-1][5] > max_end_clip or alns[i][4] > max_end_clip: continue
			key = gf.canontip(alns[i-1][3], alns[i][2])
			if key not in allowed_gaps: continue
			fw = (key == (alns[i-1][3], alns[i][2]))
			gap_node = (">" if fw else "<") + gap_names[key]
			used_gap_alns.add(alns[i-1])
			used_gap_alns.add(alns[i])
			if last_used_index is None:
				current_read_start = alns[i-1][0]
				current_read_end = alns[i][1]
				current_path = alns[i-1][7] + (gap_node,) + alns[i][7]
				current_left_clip = alns[i-1][4]
				current_right_clip = alns[i][5]
				last_used_index = i
				continue
			assert last_used_index is not None
			if last_used_index == i-1:
				assert alns[i][1] > current_read_end
				current_read_end = alns[i][1]
				current_path += (gap_node,) + alns[i][7]
				current_right_clip = alns[i][5]
				last_used_index = i
				continue
			assert last_used_index != i-1
			path_length = 0
			assert len(current_path) >= 1
			path_length = len(node_seqs[current_path[0][1:]])
			for j in range(1, len(current_path)):
				path_length += len(node_seqs[current_path[j][1:]])
				path_length -= edge_overlaps[gf.canontip(current_path[j-1], gf.revnode(current_path[j]))]
			assert path_length > current_left_clip + current_right_clip
			f.write(name + "\t" + str(alns[i][6]) + "\t" + str(current_read_start) + "\t" + str(current_read_end) + "\t+\t" + "".join(current_path) + "\t" + str(path_length) + "\t" + str(current_left_clip) + "\t" + str(path_length - current_right_clip) + "\t0\t0\t60\n")
			last_used_index = i
			current_read_start = alns[i-1][0]
			current_read_end = alns[i][1]
			current_path = alns[i-1][7] + (gap_node,) + alns[i][7]
			current_left_clip = alns[i-1][4]
			current_right_clip = alns[i][5]
		if last_used_index is not None:
			path_length = 0
			assert len(current_path) >= 1
			path_length = len(node_seqs[current_path[0][1:]])
			for j in range(1, len(current_path)):
				path_length += len(node_seqs[current_path[j][1:]])
				path_length -= edge_overlaps[gf.canontip(current_path[j-1], gf.revnode(current_path[j]))]
			if path_length > current_left_clip + current_right_clip:
				f.write(name + "\t" + str(alns[i][6]) + "\t" + str(current_read_start) + "\t" + str(current_read_end) + "\t+\t" + "".join(current_path) + "\t" + str(path_length) + "\t" + str(current_left_clip) + "\t" + str(path_length - current_right_clip) + "\t0\t0\t60\n")

with open(nongap_aln_file_out, "w") as f2:
	with open(input_aln_file) as f:
		for l in f:
			parts = l.strip().split('\t')
			readname = parts[0].split(' ')[0]
			readlen = int(parts[1])
			readstart = int(parts[2])
			readend = int(parts[3])
			path = tuple(parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t'))
			alnstart = gf.revnode(path[0])
			alnend = path[-1]
			leftclip = int(parts[7])
			rightclip = int(parts[6]) - int(parts[8])
			assert rightclip >= 0
			alnkey = (readstart, readend, alnstart, alnend, leftclip, rightclip, readlen, path)
			if alnkey not in used_gap_alns: f2.write(l)

sys.stderr.write("inserted " + str(next_gap_id-1) + " gaps\n")
