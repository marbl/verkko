#!/usr/bin/env python

import random
import sys
import re
import graph_functions as gf

mapping_file = sys.argv[1]
edge_overlap_file = sys.argv[2]

#Read alignments to the graph. Only hifi alignments to mbg and gaps are currently used
read_alignment_file = sys.argv[3]

#Either paths from rukki or just single nodes
paths_file = sys.argv[4]

nodelens_file = sys.argv[5]
layout_output = sys.argv[6]
layscf_output = sys.argv[7]

MAX_GAP_SIZE = 100000
min_contig_no_trim = 500000
min_read_len_fraction = 0.5
min_read_fromend_fraction = min_read_len_fraction/1.5
min_exact_len_fraction = min_read_len_fraction/3

#transform paths to base elements - mbg nodes and gaps.
def get_leafs(path, mapping, edge_overlaps, raw_node_lens):
	path_len = 0
	for i in range(0, len(path)):
		path_len += raw_node_lens[path[i][1:]]
		if i > 0: path_len -= edge_overlaps[gf.canon(path[i-1], path[i])]
	result = [(n, 0, raw_node_lens[n[1:]]) for n in path]
	overlaps = []
	for i in range(1, len(path)):
		overlaps.append(edge_overlaps[gf.canon(path[i-1], path[i])])
	current_len = 0
	for i in range(0, len(result)):
		assert result[i][2] > result[i][1]
		assert result[i][2] <= raw_node_lens[result[i][0][1:]]
		assert result[i][1] >= 0
		current_len += result[i][2] - result[i][1]
		if i > 0: current_len -= overlaps[i-1]
	assert current_len == path_len
	while True:
		any_replaced = False
		new_result = []
		new_overlaps = []
		for i in range(0, len(result)):
			if result[i][0][1:] not in mapping:
				new_result.append(result[i])
				if i > 0: new_overlaps.append(overlaps[i-1])
			else:
				any_replaced = True
				part = [(n, 0, raw_node_lens[n[1:]]) for n in mapping[result[i][0][1:]][0]]
				part[0] = (part[0][0], part[0][1] + mapping[result[i][0][1:]][1], part[0][2])
				part[-1] = (part[-1][0], part[-1][1], part[-1][2] - mapping[result[i][0][1:]][2])
				if result[i][0][0] == "<":
					part = [(gf.revnode(n[0]), raw_node_lens[n[0][1:]] - n[2], raw_node_lens[n[0][1:]] - n[1]) for n in part[::-1]]
				old_start_clip = result[i][1]
				old_end_clip = (raw_node_lens[result[i][0][1:]] - result[i][2])
				part[0] = (part[0][0], part[0][1] + old_start_clip, part[0][2])
				part[-1] = (part[-1][0], part[-1][1], part[-1][2] - old_end_clip)
				new_result += part
				if i > 0: new_overlaps.append(overlaps[i-1])
				for j in range(1, len(part)):
					new_overlaps.append(edge_overlaps[gf.canon(part[j-1][0], part[j][0])])
		assert len(new_result) == len(new_overlaps)+1
		assert len(new_result) >= len(result)
		if not any_replaced: break
		result = new_result
		overlaps = new_overlaps
		current_len = 0
		for i in range(0, len(result)):
			# strangely, this assertion is not always true.
			# The ONT based k-mer increase can create a node where the overlap is greater than the initial MBG node size
			# and in that case the initial MBG node will have a "negative" length within the contig
			# assert result[i][2] > result[i][1]
			assert result[i][2] <= raw_node_lens[result[i][0][1:]]
			assert result[i][1] >= 0
			current_len += result[i][2] - result[i][1]
			if i > 0: current_len -= overlaps[i-1]
		assert current_len == path_len
	return (result, overlaps)

#this gives us matches of individual read to contigs(= rukki paths or isolated utig4 node)
#path is an alignment of read (for most cases, hifi read to utig1 graph)
#node_poses - map from nodes to the list of contigs and their positions in contigs
def get_matches(path, node_poses, contig_nodeseqs, raw_node_lens, edge_overlaps, pathleftclip, pathrightclip, readleftclip, readrightclip, readlen, readstart, readend, gap):
	longest = None
	result = []
	read_path_start_pos = []
	read_path_end_pos = []
	read_node_start_pos = []
	read_node_end_pos = []
	path_start_pos = 0
	read_start_along_path = pathleftclip
	read_end_along_path = pathleftclip + readend - readstart
	for i in range(0, len(path)):
		if i > 0:
			path_start_pos += raw_node_lens[path[i-1][1:]]
			path_start_pos -= edge_overlaps[gf.canon(path[i-1], path[i])]
		path_end_pos = path_start_pos + raw_node_lens[path[i][1:]]
		read_start_pos = 0
		if path_start_pos < read_start_along_path:
			read_start_pos = readstart
			read_node_start = read_start_along_path - path_start_pos
		else:
			read_start_pos = readstart + path_start_pos - read_start_along_path
			read_node_start = 0
		read_end_pos = 0
		if path_end_pos > read_end_along_path:
			read_end_pos = readend
			read_node_end = raw_node_lens[path[i][1:]] - (path_end_pos - read_end_along_path)
		else:
			read_end_pos = readend - (read_end_along_path - path_end_pos)
			read_node_end = raw_node_lens[path[i][1:]]
		assert read_start_pos >= readstart
		if read_end_pos <= read_start_pos: return []
		assert read_end_pos > read_start_pos
		assert read_end_pos <= readend
		read_path_start_pos.append(read_start_pos)
		read_path_end_pos.append(read_end_pos)
		read_node_start_pos.append(read_node_start)
		read_node_end_pos.append(read_node_end)
	for i in range(0, len(path)):
		if path[i][1:] not in node_poses: continue
		for startpos in node_poses[path[i][1:]]:
			(contig, index, fw) = startpos
			if path[i][0] == '<': fw = not fw
			assert (not fw) or contig_nodeseqs[contig][index][0] == path[i]
			assert (fw) or gf.revnode(contig_nodeseqs[contig][index][0]) == path[i]
			node_min_start = contig_nodeseqs[contig][index][1]
			node_max_end = contig_nodeseqs[contig][index][2]
			if fw:
				if node_min_start >= read_node_end_pos[i] or node_max_end <= read_node_start_pos[i]: continue
				extra_left_clip = 0
				extra_right_clip = 0
				if node_min_start > read_node_start_pos[i]:
					extra_left_clip = node_min_start - read_node_start_pos[i]
					node_start_offset = 0
				else:
					node_start_offset = read_node_start_pos[i] - node_min_start
				if node_max_end < read_node_end_pos[i]:
					extra_right_clip = read_node_end_pos[i] - node_max_end
					node_end_offset = node_max_end
				else:
					node_end_offset = read_node_end_pos[i]
				read_start_offset = read_path_start_pos[i] + extra_left_clip
				read_end_offset = read_path_end_pos[i] - extra_right_clip
			else:
				read_wanted_start_pos = raw_node_lens[path[i][1:]] - read_node_end_pos[i]
				read_wanted_end_pos = raw_node_lens[path[i][1:]] - read_node_start_pos[i]
				if node_min_start >= read_wanted_end_pos or node_max_end <= read_wanted_start_pos: continue
				extra_left_clip = 0
				extra_right_clip = 0
				if node_min_start > read_wanted_start_pos:
					extra_right_clip = node_min_start - read_wanted_start_pos
					node_end_offset = node_max_end
				else:
					node_end_offset = node_max_end - (read_wanted_start_pos - node_min_start)
				if node_max_end < read_wanted_end_pos:
					extra_left_clip = read_wanted_end_pos - node_max_end
					node_start_offset = 0
				else:
					node_start_offset = node_max_end - read_wanted_end_pos
				read_start_offset = read_path_start_pos[i] + extra_left_clip
				read_end_offset = read_path_end_pos[i] - extra_right_clip
			if read_end_offset <= read_start_offset: continue
			assert read_start_offset >= 0
			assert read_end_offset > read_start_offset
			assert read_end_offset <= readlen
			assert extra_left_clip + extra_right_clip < read_path_end_pos[i] - read_path_start_pos[i]
			if node_end_offset <= node_start_offset: continue
			assert node_start_offset >= 0
			assert node_end_offset > node_start_offset
			assert node_end_offset <= raw_node_lens[path[i][1:]]
			match_bp_len = read_path_end_pos[i] - read_path_start_pos[i] - (extra_left_clip + extra_right_clip)
			result.append((match_bp_len, contig, index, fw, i, node_start_offset, node_end_offset, read_start_offset, read_end_offset, readlen, gap))
	return result

def get_exact_match_length(clusters):
	clusters.sort(key=lambda x: x[0])
	result = 0
	max_match = 0
	for cluster in clusters:
		start = cluster[0]
		end = cluster[1]
		if end <= max_match: continue
		if start > max_match:
			result += end - start
			max_match = end
		else:
			result += end - max_match
			max_match = end
	return result

raw_node_lens = {}
with open(nodelens_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in raw_node_lens or raw_node_lens[parts[0]] == int(parts[1])
		raw_node_lens[parts[0]] = int(parts[1])

edge_overlaps = {}
with open(edge_overlap_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] == "L"
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		overlap = int(parts[5][:-1])
		key = gf.canon(fromnode, tonode)
		if key in edge_overlaps: assert edge_overlaps[key] == overlap
		edge_overlaps[key] = overlap

node_mapping = {}
with open(mapping_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in node_mapping
		path = parts[1].split(':')[0].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		left_clip = int(parts[1].split(':')[1])
		right_clip = int(parts[1].split(':')[2])
		node_mapping[parts[0]] = (path, left_clip, right_clip)
		left_len = raw_node_lens[parts[0]]
		right_len = 0
		for i in range(0, len(path)):
			right_len += raw_node_lens[path[i][1:]]
			if i > 0: right_len -= edge_overlaps[gf.canon(path[i-1], path[i])]
		assert left_len == right_len - left_clip - right_clip

pieceid = 0

#all these contains info about contigs - here nodes or rukki paths splitted by N
#paths are transformed into mbg nodes and gaps with get_leafs procedure
contig_lens = {}
contig_nodeseqs = {}
contig_nodeoverlaps = {}
contig_node_offsets = {}
contig_pieces = {}
with open(paths_file) as f:
	for l in f:
		lp	= l.strip().split('\t')

		#  Find all words that
		#	begin with [<>], contain anything but [
		#	begin with [N, contain digits and end with N] or N:optional-description]
		#	we dump the description here and anly keep the N, digits N] part
		#
		fullname = lp[0]
		pathfull = re.findall(r"([<>][^[]+|\[N\d+N(?:[^\]]+){0,1}\])", lp[1])

		contig_pieces[fullname] = []

		for pp in pathfull:
			#pp is either path without gaps or gap. In latest case do nothing
			gp = re.match(r"\[(N\d+N)(?:[^\]]+){0,1}\]", pp)
			if gp:
				tuned_numn = min(round(int(gp.group(1)[1:-1]) * 1.5), MAX_GAP_SIZE)
				contig_pieces[fullname].append("[N" + str(tuned_numn) + "N]")
				continue

			pieceid = pieceid + 1
			pathname = f"piece{pieceid:06d}"

			(path, overlaps) = get_leafs(re.findall(r"[<>][^<>]+", pp), node_mapping, edge_overlaps, raw_node_lens)
			# skip a path if the only thing in it is a gapfill
			if len(path) == 1 and path[0][0][1:4] == "gap":
				continue
			contig_pieces[fullname].append(pathname)

			contig_nodeseqs[pathname] = path
			contig_nodeoverlaps[pathname] = overlaps
			contig_node_offsets[pathname] = []
			pos = 0
			for i in range(0, len(path)-1):
				contig_node_offsets[pathname].append(pos)
				pos += path[i][2] - path[i][1]
				pos -= overlaps[i]
			contig_node_offsets[pathname].append(pos)
			contig_lens[pathname] = contig_node_offsets[pathname][-1] + path[-1][2] - path[-1][1]
			check_len = 0
			for i in range(0, len(path)):
				check_len += path[i][2] - path[i][1]
				if i > 0: check_len -= overlaps[i-1]
			assert contig_lens[pathname] == check_len
			pathstr = ""
			for i in range(0, len(path)):
				pathstr += path[i][0] + ":" + str(path[i][1]) + ":" + str(path[i][2]) + "(" + str(contig_node_offsets[pathname][i]) + ")"
				if i < len(path)-1: pathstr += "-" + str(overlaps[i])
			# sys.stderr.write(pathname + "\t" + "".join(str(n[0]) + ":" + str(n[1]) + ":" + str(n[2]) for n in path) + "\n")
			sys.stderr.write(pathname + "\t" + pathstr + "\n")

		contig_pieces[fullname].append("end")

#map from node to list of contigs where it occures
node_poses = {}
for contigname in contig_nodeseqs:
	for i in range(0, len(contig_nodeseqs[contigname])):
		nodename = contig_nodeseqs[contigname][i][0][1:]
		if contig_nodeseqs[contigname][i][0][0] == ">":
			if nodename not in node_poses: node_poses[nodename] = []
			node_poses[nodename].append((contigname, i, True))
		else:
			if nodename not in node_poses: node_poses[nodename] = []
			node_poses[nodename].append((contigname, i, False))

read_name_to_id = {}
next_read_id = 0

matches_per_read = {}
readname_to_paths = {}
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
		gap = False
		for node in path:
			if node[1:4] == "gap":
				gap = True
				break
		readname_to_paths[readname] = [path, gap]

		matches = get_matches(path, node_poses, contig_nodeseqs, raw_node_lens, edge_overlaps, pathleftclip, pathrightclip, readleftclip, readrightclip, readlen, readstart, readend, gap)
		if len(matches) == 0: continue
		if readname not in matches_per_read: matches_per_read[readname] = []
		matches_per_read[readname] += matches

contig_contains_reads = {}
for readname in matches_per_read:
	for match in matches_per_read[readname]:
		(match_bp_size, contig, contigstart, fw, pathstart, node_start_offset, node_end_offset, readstart, readend, readlen, gap) = match
		assert readstart < readend

		if fw:
			contigpos = contig_node_offsets[contig][contigstart]
			contigpos += node_start_offset
			contigpos -= readstart

			if contig not in contig_contains_reads: contig_contains_reads[contig] = {}
			if readname not in contig_contains_reads[contig]: contig_contains_reads[contig][readname] = []
			len_readstart = readstart
			len_readend = readend
			if contigstart == 0 and node_start_offset <= 50: len_readstart = 0
			if contigstart == len(contig_node_offsets[contig]) - 1 and node_end_offset >= contig_nodeseqs[contig][contigstart][2]-50: len_readend = readlen
			if gap:
				len_readstart = 0
				len_readend = readlen
				readstart = 0
				readend = readlen
			contig_contains_reads[contig][readname].append((contigpos, contigpos + readlen, len_readstart, len_readend, readlen, readstart, readend))
		else:
			contigpos = contig_node_offsets[contig][contigstart]
			contigpos += contig_nodeseqs[contig][contigstart][2] - contig_nodeseqs[contig][contigstart][1]
			contigpos -= node_start_offset
			contigpos += readstart
			contigpos -= readlen

			if contig not in contig_contains_reads: contig_contains_reads[contig] = {}
			if readname not in contig_contains_reads[contig]: contig_contains_reads[contig][readname] = []
			len_readstart = readstart
			len_readend = readend
			if contigstart == 0 and node_end_offset >= contig_nodeseqs[contig][contigstart][2]-50: len_readend = readlen
			if contigstart == len(contig_node_offsets[contig]) - 1 and node_start_offset <= 50: len_readstart = 0
			if gap:
				len_readstart = 0
				len_readend = readlen
				readstart = 0
				readend = readlen
			contig_contains_reads[contig][readname].append((contigpos + readlen, contigpos, len_readstart, len_readend, readlen, readstart, readend))

#here we clusterize separate matches to the nodes in the path to clusters (by match position in contig)

read_clusters = {}
total_banned = 0
total_notbanned = 0
for contig in contig_contains_reads:
#	print (f"debugggging  {contig}")
	cont_seq = ""
	for i in range(0, len(contig_nodeseqs[contig])):
		cont_seq += (contig_nodeseqs[contig][i][0])
#	print (cont_seq)

	sys.stdout.flush()
#lets try to ban the reads that have a suffix or prefix that contradict to the contig
	contig_ends = [contig_nodeseqs[contig][0][0], contig_nodeseqs[contig][-1][0]]

	contig_nodes_no_dir = set()
	for i in range(0, len(contig_nodeseqs[contig])):
		contig_nodes_no_dir.add(contig_nodeseqs[contig][i][0][1:])

	for readname in contig_contains_reads[contig]:
		if readname not in read_clusters: read_clusters[readname] = []
		lines = contig_contains_reads[contig][readname]
		assert len(lines) > 0
		readlen = lines[0][4]

		#lets not ban gap containing reads
		if not readname_to_paths[readname][1]:
			to_ban = False
			near_contig_end = False
			#if read_end not in contig, then something likely went wrong..
			#can be replaced by more accurate check to exclude case where read end is somewhere else, but it should not make big sense.

			to_ban_left = False
			to_ban_right = False
			if not (readname_to_paths[readname][0][0][1:] in contig_nodes_no_dir):
				to_ban_left = True
			if not (readname_to_paths[readname][0][-1][1:] in contig_nodes_no_dir):
				to_ban_right = True
			#lets always allow read to extend contig
			contains_contig_ends = [False, False]
			for node in readname_to_paths[readname][0]:
				for i in range(0,2):
					#nodes same without respect to the orientation
					if node[1:] == contig_ends[i][1:]:
						contains_contig_ends[i] = True
						#read going different direction with contig, so should switch ban logic
						#rc loop may confuse this logic, but it would be weird anyway.
						'''if node[0] != contig_ends[i][0]:
							tmp = to_ban_left
							to_ban_left = to_ban_right
							to_ban_right = tmp
							print ("switching")
							print (readname_to_paths[readname][0]) '''


			if to_ban_left and not(contains_contig_ends[0]):
				to_ban = True
			if to_ban_right and not(contains_contig_ends[1]):
				to_ban = True
			if to_ban:
#				print(f" {to_ban_left} {to_ban_right} banning")
#				print (readname_to_paths[readname][0])
				total_banned +=1
				continue
			else:
				total_notbanned +=1

		lines.sort(key=lambda x: min(x[0], x[1]))
		fwcluster = None
		bwcluster = None
		from_end  = 0

		for line in lines:
			fw = line[1] > line[0]
			contigstart = min(line[0], line[1])
			contigend = max(line[0], line[1])
			readstart = line[2]
			readend = line[3]
			real_readstart = line[5]
			real_readend = line[6]
			assert readstart < readend
			if fw:
				if fwcluster is None:
					fwcluster = (contigstart, contigend, readstart, readend, [(real_readstart, real_readend)])
				#this is the place where alignment to different nodes are actually united
				elif contigstart < fwcluster[0] + 50 and contigend < fwcluster[1] + 50:
					fwcluster = (contigstart, contigend, min(fwcluster[2], readstart), max(fwcluster[3], readend), fwcluster[4] + [(real_readstart, real_readend)])
				else:
					#sys.stderr.write("Checking read %s in contig %s with len %s and span %s and exact %s matches fw norm %s\n"%(readname, contig, readlen, (fwcluster[3] - fwcluster[2]), get_exact_match_length(fwcluster[4]), fwcluster[4]))
					if fwcluster[3] - fwcluster[2] >= readlen * min_read_len_fraction and get_exact_match_length(fwcluster[4]) >= readlen * min_exact_len_fraction:
						read_clusters[readname].append((contig, fwcluster[0], fwcluster[1], get_exact_match_length(fwcluster[4])))
					fwcluster = (contigstart, contigend, readstart, readend, [(real_readstart, real_readend)])
			else:
				if bwcluster is None:
					bwcluster = (contigstart, contigend, readstart, readend, [(real_readstart, real_readend)])
				elif contigstart < bwcluster[0] + 50 and contigend < bwcluster[1] + 50:
					bwcluster = (contigstart, contigend, min(bwcluster[2], readstart), max(bwcluster[3], readend), bwcluster[4] + [(real_readstart, real_readend)])
				else:
					#sys.stderr.write("Checking read %s in contig %s with len %s and span %s and exact %s matches bw norm %s\n"%(readname, contig, readlen, (bwcluster[3] - bwcluster[2]), get_exact_match_length(bwcluster[4]), bwcluster[4]))
					if bwcluster[3] - bwcluster[2] >= readlen * min_read_len_fraction and get_exact_match_length(bwcluster[4]) >= readlen * min_exact_len_fraction:
						read_clusters[readname].append((contig, bwcluster[1], bwcluster[0], get_exact_match_length(bwcluster[4])))
					bwcluster = (contigstart, contigend, readstart, readend, [(real_readstart, real_readend)])
		if fwcluster is not None:
			if fwcluster[0] < 0: from_end = fwcluster[1] - (fwcluster[3] - fwcluster[2])
			if fwcluster[1] > contig_lens[contig]: from_end = contig_lens[contig] - (fwcluster[0] + (fwcluster[3] - fwcluster[2]))
			if from_end > min_read_fromend_fraction * readlen:
				total_banned += 1

			#sys.stderr.write("Checking read %s in contig %s with len %s and span %s and exact %s matches fw end %s\n"%(readname, contig, readlen, (fwcluster[3] - fwcluster[2]), get_exact_match_length(fwcluster[4]), fwcluster[4]))

			if fwcluster[3] - fwcluster[2] >= readlen * min_read_len_fraction and from_end <= min_read_fromend_fraction * readlen and get_exact_match_length(fwcluster[4]) >= readlen * min_exact_len_fraction:
				read_clusters[readname].append((contig, fwcluster[0], fwcluster[1], get_exact_match_length(fwcluster[4])))
		if bwcluster is not None:
			if bwcluster[0] < 0: from_end = bwcluster[1] - (bwcluster[3] - bwcluster[2])
			if bwcluster[1] > contig_lens[contig]: from_end = contig_lens[contig] - (bwcluster[0] + (bwcluster[3] - bwcluster[2]))
			if from_end > min_read_fromend_fraction * readlen:
				total_banned += 1

			#sys.stderr.write("Checking read %s in contig %s with len %s and span %s and exact %s matches bw end %s\n"%(readname, contig, readlen, (bwcluster[3] - bwcluster[2]), get_exact_match_length(bwcluster[4]), bwcluster[4]))
			if bwcluster[3] - bwcluster[2] >= readlen * min_read_len_fraction and from_end <= min_read_fromend_fraction * readlen and get_exact_match_length(bwcluster[4]) >= readlen * min_exact_len_fraction:
				read_clusters[readname].append((contig, bwcluster[1], bwcluster[0], get_exact_match_length(bwcluster[4])))
print (f"Banned {total_banned} allowed {total_notbanned} read-paths")
contig_actual_lines = {}
read_actual_counts = {}
read_actual_contigs = {}
#selecting best alignment to the contig, if multiple best use all of them.
for readname in read_clusters:
	longest = []
	for line in read_clusters[readname]:
		if len(longest) == 0:
			longest.append(line)
		elif line[3] > longest[0][3]:
			longest = []
			longest.append(line)
		elif line[3] == longest[0][3]:
			longest.append(line)
	for line in longest:
		if line[0] not in contig_actual_lines: contig_actual_lines[line[0]] = []
		contig_actual_lines[line[0]].append((readname, line[1], line[2]))
		read_actual_counts[readname] = read_actual_counts.get(readname, 0) + 1
		if readname not in read_actual_contigs: read_actual_contigs[readname] = set()
		read_actual_contigs[readname].add(line[0])

tig_layout_file = open(f"{layout_output}", mode="w")
scf_layout_file = open(f"{layscf_output}", mode="w")
nul_layout_file = open(f"{layscf_output}.dropped", mode="w")

#  Emit a scaffold-map entry if the contig has pieces with reads assigned.
#  Also renames contigs from the rukki name ("mat_from_{verkko_name}") to
#  "hap1-{id}", arbitrarily assigning mat to hap1.
#
#  If the contig has pieces both with and without reads assigned, emit
#  a warning (that probably nobody will see).  On the other hand, if the
#  contig actually has pieces, output the scaffold map.  (The header line
#  output from rukki looks like a contig with no pieces.)
#
no_trim = set()
nameid = 1
for contig in sorted(contig_pieces.keys()):
	npieces = ngaps = nempty = 0
	ngaps   = 0
	nempty  = 0

	for i in range(len(contig_pieces[contig])):
		line = contig_pieces[contig][i]
		if re.match(r"\[N\d+N\]", line):
			ngaps   += 1   #  Actual gap.
		elif line in contig_actual_lines:
			npieces += 1   #  Piece with reads assigned.
		elif line != "end":
			nempty  += 1   #  Piece with no reads assigned.

			#This is not normal but still may happen, switching off for now.

			#if contig_lens[line] > 100000:
			#	sys.stderr.write(f"ERROR: empty entry for piece {line} in scaffold {contig} of LARGE expected length {contig_lens[line]}\n")
			#	exit(1)				
			contig_pieces[contig][i] = "[N%dN]"%(contig_lens[line])
			if npieces > 0:
				print(f"Warning: empty entry for piece {line} in scaffold {contig} of expected length {contig_lens[line]}, using Ns {contig_pieces[contig][i]} instead", file=nul_layout_file)
	#postprocessing to avoid consecutive empty pieces
	i = 0
	new_pieces = []
	previous = "NONE"
	while i < len(contig_pieces[contig]):
		line = contig_pieces[contig][i]
		i += 1
		is_gap = re.match(r"\[N\d+N\]", line)
		if is_gap:
			#just ignoring leading gaps
			if previous == "NONE":
				continue			
			elif previous == "piece":
				new_pieces.append(line)
				previous = "gap"
			else:
				#merging consecutive gaps
				previous = "gap"
				last = new_pieces.pop()
				last_int = int(last[2:-2])
				cur_int = int(line[2:-2])
				new_pieces.append("[N%dN]"%(last_int + cur_int))
		else:
			new_pieces.append(line)
			# it can actually be 'end' but who cares then
			previous = "piece"
	#removing trailing gaps
	while len(new_pieces) > 0 and re.match(r"\[N\d+N\]", new_pieces[-1]):
		new_pieces.pop()
	contig_pieces[contig] = new_pieces

	if npieces > 0:
		if nempty > 0:
			print(f"{contig} has empty pieces.  npieces {npieces} ngaps {ngaps} nempty {nempty}", file=nul_layout_file)
		#  Decode the rukki (or non-rukki) contig name into
		#  something we'll present to the user.
		#
		#  Without rukki, contig names are:
		#	utig4-0
		#
		#  With rukki, they are:
		#	ladybug-hap1_from_utig4-1736
		#	ladybug-hap2_unused_utig4-0
		#	na_unused_utig4-1
		#  (where ladybug-hap1/2 is supplied by the user)
		#
		isna	 = re.search(r"^na_unused_(.*)$", contig)
		isunused = re.search(r"^(.*)_unused_(.*)$", contig)
		isfrom   = re.search(r"^(.*)_from_(.*)$", contig)

		if   isna:
			outname = f"unassigned-{nameid:07}"
		elif isunused:
			outname = f"{isunused.group(1)}-{nameid:07}"
		elif isfrom:
			outname = f"{isfrom.group(1)}-{nameid:07}"
		else:
			outname = f"contig-{nameid:07}"

		print(f"path {outname} {contig}", file=scf_layout_file)

		for line in contig_pieces[contig]:
			if len(contig_pieces[contig]) > 2 and (line == contig_pieces[contig][0] or line == contig_pieces[contig][-2]): no_trim.add(line)
			print(line, file=scf_layout_file)

		nameid += 1
	else:
		print(f"{contig} has no reads assigned and is not output.", file=nul_layout_file)

del nameid

# we will update the count after randomly assigning reads between their equal locations
# high count reads are output in all locations to avoid introducing coverage but simple hom nodes are randomly assigned
read_output_counts = {}
for contig in sorted(contig_actual_lines.keys()):
	if len(contig_actual_lines[contig]) == 0: continue
	assert len(contig_actual_lines[contig]) > 0
	contig_actual_lines[contig].sort(key=lambda x: min(x[1], x[2]))
	start_pos = contig_actual_lines[contig][0][1]
	end_pos = contig_actual_lines[contig][0][1]
	assigned_read_count = 0

	for line in contig_actual_lines[contig]:
		# we look at randomly placed reads and pick a single location for each one
		# if we're first we'll select a location, single value between 1 and number of times it's used
		if line[0] not in read_output_counts:
			read_output_counts[line[0]] = ("", "", 0)
			#sys.stderr.write("For reach %s which occurs %d times in contigs %s we are seing it for the first time in %s and selecting"%(line[0], read_actual_counts[line[0]], str(read_actual_contigs[line[0]]), contig))
			if len(read_actual_contigs[line[0]]) <= 1:
				#sys.stderr.write(" read is used %d times in %s contigs so placing everywhere "%(read_actual_counts[line[0]], str(read_actual_contigs[line[0]])))
				read_actual_counts[line[0]] = -1 # place everywhere
			else:
				read_actual_counts[line[0]] =  random.randint(1, read_actual_counts[line[0]])
			#sys.stderr.write(" position %d\n"%(read_actual_counts[line[0]]))
		#check if we are we right instance
		read_output_counts[line[0]] = (read_output_counts[line[0]][0], read_output_counts[line[0]][1], read_output_counts[line[0]][2] + 1)
		#sys.stderr.write("The count for read %s is now %s\n"%(line[0], str(read_output_counts[line[0]])))
		if read_actual_counts[line[0]] > 0 and read_output_counts[line[0]][2] != read_actual_counts[line[0]]: continue

		# record the contig and coordinate we will use
		read_output_counts[line[0]] = (contig, line[1], read_actual_counts[line[0]])
		start_pos = min(start_pos, line[1])
		start_pos = min(start_pos, line[2])
		end_pos = max(end_pos, line[1])
		end_pos = max(end_pos, line[2])
		assigned_read_count += 1

	print(f"tig\t{contig}", file=tig_layout_file)
	print(f"len\t{end_pos - start_pos}", file=tig_layout_file)
	if end_pos-start_pos >= min_contig_no_trim or contig in no_trim:
		print(f"trm\t1", file=tig_layout_file)
	else:
		print(f"trm\t0", file=tig_layout_file)
	print(f"rds\t{assigned_read_count}", file=tig_layout_file)

	for line in contig_actual_lines[contig]:
		bgn = line[1] - start_pos
		end = line[2] - start_pos
		#sys.stderr.write("For read %s which has been selected to be in position %s we are expecting %d\n"%(line[0], str(read_output_counts[line[0]]), read_actual_counts[line[0]]))
		if read_actual_counts[line[0]] < 0 or (read_output_counts[line[0]][0] == contig and read_output_counts[line[0]][1] == line[1]):
			print(f"{line[0]}\t{bgn}\t{end}", file=tig_layout_file)
	print(f"end", file=tig_layout_file)

tig_layout_file.close()
scf_layout_file.close()
nul_layout_file.close()


