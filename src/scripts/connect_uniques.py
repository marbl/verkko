#!/usr/bin/env python

import sys
import graph_functions as gf

graph_file = sys.argv[1]
forbidden_ends_file = sys.argv[2]
resolving_paths_file = sys.argv[3]
unique_nodes_file = sys.argv[4]
# graph to stdout

#TODO: the other instance or merge do parent[left] = right but this does the opposite, due to how edges are read?
def merge(parent, left, right):
	left = gf.find(parent, left)
	right = gf.find(parent, right)
	assert parent[left] == left
	assert parent[right] == right
	parent[left] = right

uniques = set()
with open(unique_nodes_file) as f:
	for l in f:
		uniques.add(l.strip())

end_parent = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			if ">" + parts[1] not in end_parent: end_parent[">" + parts[1]] = ">" + parts[1]
			if "<" + parts[1] not in end_parent: end_parent["<" + parts[1]] = "<" + parts[1]
			if parts[1] not in uniques: merge(end_parent, ">" + parts[1], "<" + parts[1])
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromnode not in end_parent: end_parent[fromnode] = fromnode
			if tonode not in end_parent: end_parent[tonode] = tonode
			merge(end_parent, fromnode, tonode)

resolvable_ends = set()
path_covered_nodes = set()
paths = []

with open(resolving_paths_file) as f:
	for line in f:
		l = line.strip() + '>'
		last_break = 0
		path = []
		for i in range(1, len(l)):
			if l[i] == '<' or l[i] == '>':
				path.append(l[last_break:i])
				last_break = i
		if len(path) == 0: continue
		resolvable_ends.add(path[0])
		resolvable_ends.add(gf.revnode(path[-1]))
		paths.append(path)
		assert gf.find(end_parent, path[0]) == gf.find(end_parent, gf.revnode(path[-1]))

resolved_tangles = set()

with open(forbidden_ends_file) as f:
	for l in f:
		parts = l.strip().split(',')
		for part in parts:
			if part in resolvable_ends: resolvable_ends.remove(part)

for tip in resolvable_ends:
	resolved_tangles.add(gf.find(end_parent, tip))

for path in paths:
	if path[0] not in resolvable_ends or gf.revnode(path[-1]) not in resolvable_ends: continue
	for node in path[1:-1]:
		path_covered_nodes.add(node[1:])

sys.stderr.write(str(len(resolvable_ends)) + " resolvable ends\n")
sys.stderr.write(str(resolvable_ends) + "\n")
sys.stderr.write(str(len(path_covered_nodes)) + " path covered nodes\n")

node_seq = {}
base_overlaps = {}

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			node_seq[parts[1]] = parts[2]
			if parts[1] in path_covered_nodes: continue
			if parts[1] not in uniques and gf.find(end_parent, ">" + parts[1]) in resolved_tangles: continue
		elif parts[0] == 'L':
			check_from = (">" if parts[2] == "+" else "<") + parts[1]
			check_to = ("<" if parts[4] == "+" else ">") + parts[3]
			key = gf.canontip(check_from, check_to)
			base_overlaps[key] = parts[5]
			if check_from in resolvable_ends or check_to in resolvable_ends:
				continue
			if parts[1] in path_covered_nodes or parts[3] in path_covered_nodes: continue
			assert gf.find(end_parent, check_from) == gf.find(end_parent, check_to)
			if (parts[1] not in uniques or parts[3] not in uniques) and gf.find(end_parent, check_from) in resolved_tangles: continue
		print(l.strip())

# sys.stderr.write(str(len(node_seq)) + " nodes\n")

total_num_insertions = {}
with open(resolving_paths_file) as f:
	for line in f:
		path = line.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		for node in path[1:-1]:
			if node[1:] not in total_num_insertions: total_num_insertions[node[1:]] = 0
			total_num_insertions[node[1:]] += 1

num_insertions = {}

with open(resolving_paths_file) as f:
	for line in f:
		l = line.strip() + '>'
		last_break = 0
		path = []
		for i in range(1, len(l)):
			if l[i] == '<' or l[i] == '>':
				path.append(l[last_break:i])
				last_break = i
		sys.stderr.write("path " + str(path) + "\n")
		assert len(path) >= 2
		if len(path) == 0: continue
		if path[0] not in resolvable_ends or ((">" if path[-1][0] == "<" else "<") + path[-1][1:]) not in resolvable_ends: continue
		sys.stderr.write("insert path " + str(path) + "\n")
		last_node = path[0]
		for i in range(1, len(path)-1):
			if path[i][1:] not in num_insertions: num_insertions[path[i][1:]] = 0
			num_insertions[path[i][1:]] += 1
			this_node = path[i][0] + path[i][1:] + "_" + str(num_insertions[path[i][1:]])
			if total_num_insertions[path[i][1:]] == 1: this_node = path[i]
			print("S\t" + this_node[1:] + "\t" + node_seq[path[i][1:]])
			overlap = "0M"
			key = gf.canontip(path[i-1], gf.revnode(path[i]))
			assert key in base_overlaps
			overlap = base_overlaps[key]
			print("L\t" + last_node[1:] + "\t" + ("+" if last_node[0] == ">" else "-") + "\t" + this_node[1:] + "\t" + ("+" if this_node[0] == ">" else "-") + "\t" + overlap)
			last_node = this_node
		key = gf.canontip(path[-2], gf.revnode(path[-1]))
		assert key in base_overlaps
		overlap = base_overlaps[key]
		print("L\t" + last_node[1:] + "\t" + ("+" if last_node[0] == ">" else "-") + "\t" + path[-1][1:] + "\t" + ("+" if path[-1][0] == ">" else "-") + "\t" + overlap)
