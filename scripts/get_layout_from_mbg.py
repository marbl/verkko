#!/usr/bin/python

import sys

mapping_file = sys.argv[1]
edge_overlap_file = sys.argv[2]
read_alignment_file = sys.argv[3]
paths_file = sys.argv[4]
nodelens_file = sys.argv[5]
# layout to stdout

min_read_len_fraction = 0.5

def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(left, right):
	if revnode(right) + revnode(left) < left + right:
		return (revnode(right), revnode(left))
	return (left, right)

def get_leafs(path, mapping, edge_overlaps, raw_node_lens):
	path_len = 0
	for i in range(0, len(path)):
		path_len += raw_node_lens[path[i][1:]]
		if i > 0: path_len -= edge_overlaps[canon(path[i-1], path[i])]
	result = path
	overlaps = []
	for i in range(1, len(path)):
		overlaps.append(edge_overlaps[canon(path[i-1], path[i])])
	current_len = 0
	for i in range(0, len(result)):
		current_len += raw_node_lens[result[i][1:]]
		if i > 0: current_len -= overlaps[i-1]
	assert current_len == path_len
	while True:
		any_replaced = False
		new_result = []
		new_overlaps = []
		for i in range(0, len(result)):
			if result[i][1:] not in mapping:
				new_result.append(result[i])
				if i > 0: new_overlaps.append(overlaps[i-1])
			else:
				any_replaced = True
				part = mapping[result[i][1:]]
				if result[i][0] == "<":
					part = [revnode(n) for n in part[::-1]]
				new_result += part
				if i > 0: new_overlaps.append(overlaps[i-1])
				for j in range(1, len(part)):
					new_overlaps.append(edge_overlaps[canon(part[j-1], part[j])])
		assert len(new_result) == len(new_overlaps)+1
		assert len(new_result) >= len(result)
		if not any_replaced: break
		result = new_result
		overlaps = new_overlaps
		current_len = 0
		for i in range(0, len(result)):
			current_len += raw_node_lens[result[i][1:]]
			if i > 0: current_len -= overlaps[i-1]
		assert current_len == path_len
	return (result, overlaps)

def get_match_len(path, pathstart, contigpath, contigstart, fw):
	match_len = 0
	while True:
		if pathstart < 0 or pathstart == len(path): break
		if contigstart < 0 or contigstart == len(contigpath): break
		if fw and path[pathstart] != contigpath[contigstart]: break
		if not fw and path[pathstart] != revnode(contigpath[contigstart]): break
		match_len += 1
		pathstart += 1
		if fw:
			contigstart += 1
		else:
			contigstart -= 1
	return match_len

def get_path_len(path, start, end, raw_node_lens, edge_overlaps):
	result = raw_node_lens[path[start][1:]]
	for i in range(1, len(path)):
		assert raw_node_lens[path[i][1:]] > edge_overlaps[canon(path[i-1], path[i])]
		result += raw_node_lens[path[i][1:]] - edge_overlaps[canon(path[i-1], path[i])]
	return result

def get_matches(path, node_poses, contig_nodeseqs, raw_node_lens, edge_overlaps, pathleftclip, pathrightclip, readleftclip, readrightclip, readlen, readstart, readend):
	longest = None
	result = []
	for i in range(0, len(path)):
		if path[i][1:] not in node_poses: continue
		for startpos in node_poses[path[i][1:]]:
			(contig, index, fw) = startpos
			if path[i][0] == '<': fw = not fw
			assert (not fw) or contig_nodeseqs[contig][index] == path[i]
			assert (fw) or revnode(contig_nodeseqs[contig][index]) == path[i]
			match_node_len = get_match_len(path, i, contig_nodeseqs[contig], index, fw)
			assert match_node_len >= 1
			assert i + match_node_len <= len(path)
			match_bp_len = get_path_len(path, i, i + match_node_len, raw_node_lens, edge_overlaps)
			if i == 0: match_bp_len -= pathleftclip
			if i + match_node_len == len(path): match_bp_len -= pathrightclip
			result.append((match_bp_len, contig, index, fw, i, match_node_len, pathleftclip, pathrightclip, readleftclip, readrightclip, path, readlen, readstart, readend))
	return result

node_mapping = {}
with open(mapping_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in node_mapping
		path = parts[1].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		node_mapping[parts[0]] = path

edge_overlaps = {}
with open(edge_overlap_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] == "L"
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		overlap = int(parts[5][:-1])
		key = canon(fromnode, tonode)
		if key in edge_overlaps: assert edge_overlaps[key] == overlap
		edge_overlaps[key] = overlap

raw_node_lens = {}
with open(nodelens_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in raw_node_lens or raw_node_lens[parts[0]] == int(parts[1])
		raw_node_lens[parts[0]] = int(parts[1])

contig_lens = {}
contig_nodeseqs = {}
contig_nodeoverlaps = {}
contig_node_offsets = {}
with open(paths_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		pathname = parts[0]
		path = parts[1].replace('<', '\t<').replace('>', '\t>').strip().split('\t')
		(path, overlaps) = get_leafs(path, node_mapping, edge_overlaps, raw_node_lens)
		sys.stderr.write(pathname + "\t" + "".join(path) + "\n")
		contig_nodeseqs[pathname] = path
		contig_nodeoverlaps[pathname] = overlaps
		contig_node_offsets[pathname] = []
		contig_node_offsets[pathname].append(0)
		for i in range(1, len(path)):
			prev = path[i-1]
			curr = path[i]
			overlap = overlaps[i-1]
			nodelen = raw_node_lens[prev[1:]]
			contig_node_offsets[pathname].append(contig_node_offsets[pathname][-1] + nodelen - overlap)
		contig_lens[pathname] = contig_node_offsets[pathname][-1] + raw_node_lens[path[-1][1:]]

node_poses = {}
for contigname in contig_nodeseqs:
	for i in range(0, len(contig_nodeseqs[contigname])):
		nodename = contig_nodeseqs[contigname][i][1:]
		if contig_nodeseqs[contigname][i][0] == ">":
			if nodename not in node_poses: node_poses[nodename] = []
			node_poses[nodename].append((contigname, i, True))
		else:
			if nodename not in node_poses: node_poses[nodename] = []
			node_poses[nodename].append((contigname, i, False))

read_name_to_id = {}
next_read_id = 0

matches_per_read = {}
with open(read_alignment_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		readname = parts[0].split(' ')[0]
		if readname not in read_name_to_id:
			read_name_to_id[readname] = next_read_id
			next_read_id += 1
		readlen = int(parts[1])
		readleftclip = int(parts[2])
		readrightclip = int(parts[1]) - int(parts[3])
		readstart = int(parts[2])
		readend = int(parts[3])
		if not readstart < readend: print(l)
		assert readstart < readend
		pathleftclip = int(parts[7])
		pathrightclip = int(parts[6]) - int(parts[8])
		path = parts[5].replace('>', "\t>").replace('<', "\t<").strip().split('\t')
		matches = get_matches(path, node_poses, contig_nodeseqs, raw_node_lens, edge_overlaps, pathleftclip, pathrightclip, readleftclip, readrightclip, readlen, readstart, readend)
		if len(matches) == 0: continue
		if readname not in matches_per_read: matches_per_read[readname] = []
		matches_per_read[readname] += matches

contig_contains_reads = {}
for readname in matches_per_read:
	for match in matches_per_read[readname]:
		(match_bp_size, contig, contigstart, fw, pathstart, matchlen, pathleftclip, pathrightclip, readleftclip, readrightclip, path, readlen, readstart, readend) = match
		assert readstart < readend
		if fw:
			contigpos = contig_node_offsets[contig][contigstart]
			for i in range(0, pathstart):
				contigpos -= raw_node_lens[path[i][1:]] - edge_overlaps[canon(path[i], path[i+1])]
			contigpos += pathleftclip
			contigpos -= readleftclip
			if contig not in contig_contains_reads: contig_contains_reads[contig] = {}
			if readname not in contig_contains_reads[contig]: contig_contains_reads[contig][readname] = []
			len_readstart = readstart
			len_readend = readend
			if contigstart == 0: len_readstart = 0
			if contigstart + matchlen == len(contig_nodeseqs[contig]): len_readend = readlen
			if len(path) == 1 and path[0][1:4] == "gap":
				len_readstart = 0
				len_readend = readlen
			contig_contains_reads[contig][readname].append((contigpos, contigpos + readlen, len_readstart, len_readend, readlen, readstart, readend))
		else:
			pathstart += matchlen - 1
			contigstart -= matchlen - 1
			assert pathstart < len(path)
			assert contigstart >= 0
			contigpos = contig_node_offsets[contig][contigstart]
			for i in range(pathstart+1, len(path)):
				contigpos -= raw_node_lens[path[i][1:]] - edge_overlaps[canon(path[i-1], path[i])]
			contigpos += pathrightclip
			contigpos -= readrightclip
			if contig not in contig_contains_reads: contig_contains_reads[contig] = {}
			if readname not in contig_contains_reads[contig]: contig_contains_reads[contig][readname] = []
			len_readstart = readstart
			len_readend = readend
			if contigstart == len(contig_nodeseqs[contig]) - 1: len_readstart = 0
			if contigstart - matchlen == 0: len_readend = readlen
			if len(path) == 1 and path[0][1:4] == "gap":
				len_readstart = 0
				len_readend = readlen
			contig_contains_reads[contig][readname].append((contigpos + readlen, contigpos, len_readstart, len_readend, readlen, readstart, readend))

read_clusters = {}
for contig in contig_contains_reads:
	for readname in contig_contains_reads[contig]:
		if readname not in read_clusters: read_clusters[readname] = []
		lines = contig_contains_reads[contig][readname]
		assert len(lines) > 0
		readlen = lines[0][4]
		lines.sort(key=lambda x: min(x[0], x[1]))
		fwcluster = None
		bwcluster = None
		for line in lines:
			fw = line[1] > line[0]
			contigstart = min(line[0], line[1])
			contigend = max(line[0], line[1])
			readstart = line[2]
			readend = line[3]
			real_readstart = line[4]
			real_readend = line[5]
			assert readstart < readend
			if fw:
				if fwcluster is None:
					fwcluster = (contigstart, contigend, readstart, readend, real_readstart, real_readend)
				elif contigstart < fwcluster[0] + 50 and contigend < fwcluster[1] + 50:
					fwcluster = (contigstart, contigend, min(fwcluster[2], readstart), max(fwcluster[3], readend), min(fwcluster[4], real_readstart), max(fwcluster[5], real_readend))
				else:
					if fwcluster[3] - fwcluster[2] >= readlen * min_read_len_fraction: read_clusters[readname].append((contig, fwcluster[0], fwcluster[1], fwcluster[4], fwcluster[5]))
					fwcluster = (contigstart, contigend, readstart, readend, real_readstart, real_readend)
			else:
				if bwcluster is None:
					bwcluster = (contigstart, contigend, readstart, readend, real_readstart, real_readend)
				elif contigstart < bwcluster[0] + 50 and contigend < bwcluster[1] + 50:
					bwcluster = (contigstart, contigend, min(bwcluster[2], readstart), max(bwcluster[3], readend), min(bwcluster[4], real_readstart), max(bwcluster[5], real_readend))
				else:
					if bwcluster[3] - bwcluster[2] >= readlen * min_read_len_fraction: read_clusters[readname].append((contig, bwcluster[1], bwcluster[0], bwcluster[4], bwcluster[5]))
					bwcluster = (contigstart, contigend, readstart, readend, real_readstart, real_readend)
		if fwcluster is not None:
			if fwcluster[3] - fwcluster[2] >= readlen * min_read_len_fraction: read_clusters[readname].append((contig, fwcluster[0], fwcluster[1], fwcluster[4], fwcluster[5]))
		if bwcluster is not None:
			if bwcluster[3] - bwcluster[2] >= readlen * min_read_len_fraction: read_clusters[readname].append((contig, bwcluster[1], bwcluster[0], bwcluster[4], bwcluster[5]))

contig_actual_lines = {}
for readname in read_clusters:
	longest = []
	for line in read_clusters[readname]:
		if len(longest) == 0:
			longest.append(line)
		elif line[4] - line[3] > longest[0][4] - longest[0][3]:
			longest = []
			longest.append(line)
		elif line[4] - line[3] == longest[0][4] - longest[0][3]:
			longest.append(line)
	for line in longest:
		if line[0] not in contig_actual_lines: contig_actual_lines[line[0]] = []
		contig_actual_lines[line[0]].append((readname, line[1], line[2]))

for contig in contig_actual_lines:
	if len(contig_actual_lines[contig]) == 0: continue
	assert len(contig_actual_lines[contig]) > 0
	contig_actual_lines[contig].sort(key=lambda x: min(x[1], x[2]))
	start_pos = contig_actual_lines[contig][0][1]
	end_pos = contig_actual_lines[contig][0][1]
	for line in contig_actual_lines[contig]:
		start_pos = min(start_pos, line[1])
		start_pos = min(start_pos, line[2])
		end_pos = max(end_pos, line[1])
		end_pos = max(end_pos, line[2])
	print("tig\t" + contig)
	print("len\t" + str(end_pos - start_pos))
	print("rds\t" + str(len(contig_actual_lines[contig])))
	for line in contig_actual_lines[contig]:
		print(line[0] + "\t" + str(line[1] - start_pos) + "\t" + str(line[2] - start_pos))
	print("end")
