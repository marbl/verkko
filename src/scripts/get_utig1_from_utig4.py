#!/usr/bin/env python

import random
import sys
import re
import graph_functions as gf

MAX_GAP_SIZE=100000
mapping_file = sys.argv[1]
edge_overlap_file = sys.argv[2]

#Either paths from rukki or just single nodes
paths_file = sys.argv[3]

nodelens_file = sys.argv[4]

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
cut_mapping  = {}
with open(mapping_file) as f:
	for l in f:

		parts = l.strip().split('\t')
		assert parts[0] not in node_mapping
		if not re.search(r"utig\d+[a-z]?-" , parts[1]):
			continue
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

		# save the mapping of cut nodes to their respective coordinates so we can find them later
		if (len(path) == 1 and path[0][1:] in raw_node_lens):
			new_name = path[0][1:] + ":" + str(left_clip) + ":" + str(raw_node_lens[path[0][1:]]-right_clip)
			cut_mapping[new_name] = parts[0]
pieceid = 0

#all these contains info about contigs - here nodes or rukki paths splitted by N
#paths are transformed into mbg nodes and gaps with get_leafs procedure
contig_lens = {}
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
				contig_pieces[fullname].append("[N" + str(tuned_numn) + "N:gap]")
				continue

			pieceid = pieceid + 1
			pathname = fullname

			(path, overlaps) = get_leafs(re.findall(r"[<>][^<>]+", pp), node_mapping, edge_overlaps, raw_node_lens)
			# skip a path if the only thing in it is a gapfill
			if len(path) == 1 and path[0][0][1:4] == "gap":
				continue

			contig_node_offsets[pathname] = []
			pos = 0
			end = -1
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
				# build a name using the contig without the <> but also append coordinates if it's partial match to check for cut node
				# if a cut version exists, use that name instead, otherwise use the original node name
				new_name = path[i][0][1:]
				if path[i][1] != 0 or path[i][2] != raw_node_lens[new_name]:
					if path[i][0][0] == ">":
						new_name = path[i][0][1:] + ":" + str(path[i][1]) + ":" + str(path[i][2])
					else:
						new_name = path[i][0][1:] + ":" + str(raw_node_lens[new_name]-path[i][2]) + ":" + str(raw_node_lens[new_name]-path[i][1])
					if new_name not in cut_mapping:
						new_name = path[i][0][1:]

				# when we see the name in our path already and the offset is earlier than the largest we have already seen, this is an overlap
				# we skip these overlapping nodes from the path and continue at the new unique/larger offset node
				#sys.stderr.write("Checking node %s with coordinates %d-%d and offset is %d vs %d and is already used is %d\n"%(path[i][0], path[i][1], path[i][2], (contig_node_offsets[pathname][i]-path[i][1]), end, (new_name in pathstr)))
				if (contig_node_offsets[pathname][i]-path[i][1]) <= end and new_name in pathstr:
					continue
				end = contig_node_offsets[pathname][i]-path[i][1]
				if path[i][1] != 0 or path[i][2] != raw_node_lens[path[i][0][1:]]:
					if (new_name in cut_mapping):
						pathstr += path[i][0] + "_" + cut_mapping[new_name].strip().split("_")[-1]
					else:
						pathstr += path[i][0]
				else:
					pathstr += path[i][0]
			contig_pieces[fullname].append(pathstr)

for fullname in contig_pieces:		
	if "name" in fullname: continue
	sys.stdout.write(fullname + "\t" + "".join(contig_pieces[fullname]) + "\n")
