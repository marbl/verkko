#!/usr/bin/python

import sys

graph_file = sys.argv[1]
paths_file = sys.argv[2]
old_uniques_file = sys.argv[3]
min_safe_coverage = int(sys.argv[4])
max_crosslink_coverage_fraction = float(sys.argv[5])
max_crosslink_coverage_absolute = int(sys.argv[6])
# new uniques to stdout

def find(parent, n):
	while parent[parent[n]] != parent[n]: parent[n] = parent[parent[n]]
	return parent[n]

def merge(parent, left, right):
	left = find(parent, left)
	right = find(parent, right)
	assert parent[left] == left
	assert parent[right] == right
	parent[right] = left

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def tangle_is_valid(tips, paths, path_index):
	global min_safe_coverage
	global max_crosslink_coverage
	if len(tips) == 0: return False
	if len(tips) % 2 == 1: return False
	tip_connections = {}
	relevant_paths = set()
	for tip in tips:
		if tip[1:] not in path_index: continue
		for i in path_index[tip[1:]]:
			relevant_paths.add(i)
	for i in relevant_paths:
		path = paths[i]
		last_tip = None
		for i in range(0, len(path)):
			if last_tip is not None and revnode(path[i]) in tips:
				this_tip = revnode(path[i])
				if last_tip not in tip_connections: tip_connections[last_tip] = {}
				if this_tip not in tip_connections[last_tip]: tip_connections[last_tip][this_tip] = 0
				tip_connections[last_tip][this_tip] += 1
				if this_tip not in tip_connections: tip_connections[this_tip] = {}
				if last_tip not in tip_connections[this_tip]: tip_connections[this_tip][last_tip] = 0
				tip_connections[this_tip][last_tip] += 1
			if path[i] in tips: last_tip = path[i]
	tip_unique_connection = {}
	biggest_crosslink_coverage = 0
	smallest_validlink_coverage = None
	for tip in tips:
		if tip not in tip_connections:
			return False
		max_connection = (None, 0)
		for other_tip in tip_connections[tip]:
			if tip_connections[tip][other_tip] > max_connection[1]:
				max_connection = (other_tip, tip_connections[tip][other_tip])
		if max_connection[1] < min_safe_coverage:
			return False
		if smallest_validlink_coverage is None: smallest_validlink_coverage = max_connection[1]
		smallest_validlink_coverage = min(smallest_validlink_coverage, max_connection[1])
		for other_tip in tip_connections[tip]:
			if other_tip == max_connection[0]: continue
			if tip_connections[tip][other_tip] > max_crosslink_coverage_absolute: return False
			biggest_crosslink_coverage = max(biggest_crosslink_coverage, tip_connections[tip][other_tip])
		tip_unique_connection[tip] = max_connection[0]
	if biggest_crosslink_coverage >= int(smallest_validlink_coverage * max_crosslink_coverage_fraction): return False
	for tip in tip_unique_connection:
		if tip_unique_connection[tip] not in tip_unique_connection:
			return False
		if tip_unique_connection[tip_unique_connection[tip]] != tip:
			return False
	return True

def remove_contained(tips):
	result = set()
	for tip in tips:
		if revnode(tip) not in tips: result.add(tip)
	return result

def get_invalid_tangle_group(start_tangle, tangle_parent, tangle_uniques, valid_tangles):
	assert start_tangle not in valid_tangles
	tangles = set()
	stack = []
	stack.append(start_tangle)
	while len(stack) >= 1:
		top = stack.pop()
		if top in tangles: continue
		tangles.add(top)
		for tip in tangle_uniques[top]:
			if tangle_parent[revnode(tip)] not in valid_tangles:
				stack.append(tangle_parent[revnode(tip)])
	return tuple(tangles)

old_uniques = set()
with open(old_uniques_file) as f:
	for l in f:
		old_uniques.add(l.strip())

paths = []
with open(paths_file) as f:
	for l in f:
		path = l.strip().replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		fix_path = []
		for node in path:
			if node[1:] in old_uniques: fix_path.append(node)
		if len(fix_path) >= 2: paths.append(fix_path)

path_index = {}
for i in range(0, len(paths)):
	nodes = set()
	for node in paths[i]: nodes.add(node[1:])
	for node in nodes:
		if node not in path_index: path_index[node] = []
		path_index[node].append(i)

tangle_parent = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			if ">" + parts[1] not in tangle_parent: tangle_parent[">" + parts[1]] = ">" + parts[1]
			if "<" + parts[1] not in tangle_parent: tangle_parent["<" + parts[1]] = "<" + parts[1]
			if parts[1] not in old_uniques:
				merge(tangle_parent, ">" + parts[1], "<" + parts[1])
		if parts[0] == "L":
			fromtip = (">" if parts[2] == "+" else "<") + parts[1]
			totip = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromtip not in tangle_parent: tangle_parent[fromtip] = fromtip
			if totip not in tangle_parent: tangle_parent[totip] = totip
			merge(tangle_parent, fromtip, totip)

tangle_uniques = {}
for n in tangle_parent:
	if n[1:] not in old_uniques: continue
	key = find(tangle_parent, n)
	if key not in tangle_uniques: tangle_uniques[key] = set()
	tangle_uniques[key].add(n)

removed_nodes = set()
valid_tangles = set()
for tangle in tangle_uniques:
	if tangle_is_valid(tangle_uniques[tangle], paths, path_index):
		valid_tangles.add(tangle)
		continue
	filtered_tangle = remove_contained(tangle_uniques[tangle])
	if tangle_is_valid(filtered_tangle, paths, path_index):
		for node in tangle_uniques[tangle]:
			if node not in filtered_tangle:
				removed_nodes.add(node[1:])
		valid_tangles.add(tangle)

invalid_tangle_groups = set()
for tangle in tangle_uniques:
	if tangle in valid_tangles: continue
	invalid_tangle_groups.add(get_invalid_tangle_group(tangle, tangle_parent, tangle_uniques, valid_tangles))

for invalid_tangle in invalid_tangle_groups:
	tangle_tips = set()
	for tangle in invalid_tangle:
		for tip in tangle_uniques[tangle]:
			tangle_tips.add(tip)
	tangle_boundary = remove_contained(tangle_tips)
	if tangle_is_valid(tangle_boundary, paths, path_index):
		for tip in tangle_tips:
			if tip not in tangle_boundary:
				removed_nodes.add(tip[1:])

sys.stderr.write("removed " + str(len(removed_nodes)) + " nodes\n")
for node in removed_nodes:
	sys.stderr.write(node + "\n")

for node in old_uniques:
	if node not in removed_nodes:
		print(node)
