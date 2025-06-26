#!/usr/bin/env python

import sys
import graph_functions as gf

graph_file = sys.argv[1]
mapping_file = sys.argv[2]
edge_overlap_file = sys.argv[3]
nodelens_file = sys.argv[4]
original_coverage_csv = sys.argv[5]
# output coverage csv to stdout

def get_leafs(path, mapping, edge_overlaps, raw_node_lens, coverages):
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
			if result[i][0][1:] not in mapping or result[i][0][1:] in coverages:
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

mapping = {}
with open(mapping_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		assert parts[0] not in mapping
		path = parts[1].split(':')[0].replace('<', "\t<").replace('>', "\t>").strip().split('\t')
		left_clip = int(parts[1].split(':')[1])
		right_clip = int(parts[1].split(':')[2])
		mapping[parts[0]] = (path, left_clip, right_clip)

original_coverages = {}
with open(original_coverage_csv) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node":
			assert(parts[1] == "coverage" and parts[2] == "length")
			continue
		original_coverages[parts[0]] = float(parts[1])

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

nodes_per_contig = {}
contig_node_offsets = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			(path, overlaps) = get_leafs([">" + parts[1]], mapping, edge_overlaps, raw_node_lens, original_coverages)
			contig_node_offsets[parts[1]] = []
			pos = 0
			for i in range(0, len(path)-1):
				contig_node_offsets[parts[1]].append(pos)
				pos += path[i][2] - path[i][1]
				pos -= overlaps[i]
			contig_node_offsets[parts[1]].append(pos)
			nodes_per_contig[parts[1]] = path

node_occurrences = {}
for contig in nodes_per_contig:
	for i in range(0, len(nodes_per_contig[contig])):
		node = nodes_per_contig[contig][i][0][1:]
		start_pos = nodes_per_contig[contig][i][1]
		end_pos = nodes_per_contig[contig][i][2]
		if end_pos <= start_pos: continue
		if nodes_per_contig[contig][i][0][0] == "<":
			start_pos = raw_node_lens[node] - start_pos
			end_pos = raw_node_lens[node] - end_pos
			(start_pos, end_pos) = (end_pos, start_pos)
		assert start_pos < end_pos
		if node not in node_occurrences: node_occurrences[node] = []
		node_occurrences[node].append((start_pos, end_pos, (contig, contig_node_offsets[contig][i], nodes_per_contig[contig][i][0][0])))

split_forbidden_intervals = {}
for node in node_occurrences:
	occurrences = list(set(node_occurrences[node]))
	occurrences.sort(key = lambda x: x[0])
	max_end = (occurrences[0][1], occurrences[0][2])
	second_max_end = (occurrences[0][1], occurrences[0][2])
	for i in range(1, len(occurrences)):
		if occurrences[i][0] <= max_end[0] and occurrences[i][2] != max_end[1]:
			if node not in split_forbidden_intervals: split_forbidden_intervals[node] = []
			split_forbidden_intervals[node].append((occurrences[i][0], min(max_end[0], occurrences[i][1])))
		if occurrences[i][0] <= second_max_end[0] and occurrences[i][2] != second_max_end[1]:
			if node not in split_forbidden_intervals: split_forbidden_intervals[node] = []
			split_forbidden_intervals[node].append((occurrences[i][0], min(second_max_end[0], occurrences[i][1])))
		if occurrences[i][2] == max_end[1]:
			max_end = (max(max_end[0], occurrences[i][1]), max_end[1])
		elif occurrences[i][2] == second_max_end[1]:
			second_max_end = (max(second_max_end[0], occurrences[i][1]), second_max_end[1])
		elif occurrences[i][1] > max_end[0]:
			assert second_max_end[0] <= max_end[0]
			second_max_end = max_end
			max_end = (occurrences[i][1], occurrences[i][2])
		elif occurrences[i][1] > second_max_end[0]:
			assert second_max_end[0] <= max_end[0]
			assert occurrences[i][1] <= max_end[0]
			second_max_end = (occurrences[i][1], occurrences[i][2])

forbidden_intervals = {}
for node in split_forbidden_intervals:
	split_forbidden_intervals[node].sort(key=lambda x: x[0])
	start = 0
	end = 0
	for interval in split_forbidden_intervals[node]:
		if interval[0] > end:
			if end > start:
				if node not in forbidden_intervals: forbidden_intervals[node] = []
				forbidden_intervals[node].append((start, end))
			start = interval[0]
			end = interval[1]
		else:
			end = max(end, interval[1])
	if end > start:
		if node not in forbidden_intervals: forbidden_intervals[node] = []
		forbidden_intervals[node].append((start, end))

print("node\tcoverage\tlength")
for contig in nodes_per_contig:
	length_sum = 0
	coverage_sum = 0
	for triple in nodes_per_contig[contig]:
		node = triple[0][1:]
		start = triple[1]
		end = triple[2]
		if triple[0][0] == "<":
			start = raw_node_lens[node] - start
			end = raw_node_lens[node] - end
			(start, end) = (end, start)
		if node not in original_coverages: continue
		if node in forbidden_intervals:
			for interval in forbidden_intervals[node]:
				assert interval[1] > interval[0]
				assert interval[0] <= start or interval[1] >= end
				if interval[0] > start:
					end = min(end, interval[0])
				else:
					start = max(start, interval[1])
		if start >= end: continue
		length_sum += end - start
		coverage_sum += original_coverages[node] * (end - start)
	if length_sum == 0: continue
	coverage = float(coverage_sum) / float(length_sum)
	print(contig + "\t" + "{:.2f}".format(coverage) + "\t%d"%(length_sum))
