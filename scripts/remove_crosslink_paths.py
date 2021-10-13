#!/usr/bin/python

import sys

unique_nodes_file = sys.argv[1]
picked_paths_file = sys.argv[2]
counted_paths_file = sys.argv[3]
# paths to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:]

def startnode(p):
	path = p.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	return path[0]

def endnode(p):
	path = p.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	return revnode(path[-1])

def check_side(start):
	global forbidden_paths
	if start not in paths_per_unique: return
	if len(paths_per_unique[start]) != 2: return
	other_side = set()
	used_paths = set()
	for path in paths_per_unique[start]:
		other_side.add(path[1])
		used_paths.add(path[0])
	if len(other_side) != 2: return
	has_two = False
	this_side = set()
	has_many = False
	for side in other_side:
		for path in paths_per_unique[side]:
			if len(paths_per_unique[side]) == 2: has_two = True
			if len(paths_per_unique[side]) > 2: has_many = True
			for side2 in paths_per_unique[side]:
				this_side.add(side2[1])
				used_paths.add(side2[0])
				if len(paths_per_unique[side2[1]]) > 2: has_many = True
	if not has_two: return
	if has_many: return
	if len(this_side) != 2: return
	total_used_nodes = set()
	for node in this_side: total_used_nodes.add(node[1:])
	for node in other_side: total_used_nodes.add(node[1:])
	if len(total_used_nodes) != 4: return
	if len(used_paths) < 3: return
	if len(used_paths) > 4: return
	abundant_paths = list(used_paths)
	abundant_paths.sort(key=lambda x: -path_counts[x])
	assert path_counts[abundant_paths[1]] >= path_counts[abundant_paths[2]]
	if path_counts[abundant_paths[1]] == path_counts[abundant_paths[2]]: return
	used_here = set()
	used_here.add(startnode(abundant_paths[0]))
	used_here.add(startnode(abundant_paths[0]))
	used_here.add(endnode(abundant_paths[1]))
	used_here.add(endnode(abundant_paths[1]))
	for i in range(2, len(abundant_paths)): forbidden_paths.add(abundant_paths[i])


unique_nodes = set()
with open(unique_nodes_file) as f:
	for l in f:
		unique_nodes.add(l.strip())

paths_per_unique = {}
all_paths = set()
with open(picked_paths_file) as f:
	for l in f:
		path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		all_paths.add(l.strip())
		if path[0] not in paths_per_unique: paths_per_unique[path[0]] = set()
		paths_per_unique[path[0]].add((l.strip(), revnode(path[-1])))
		if revnode(path[-1]) not in paths_per_unique: paths_per_unique[revnode(path[-1])] = set()
		paths_per_unique[revnode(path[-1])].add((l.strip(), path[0]))

path_counts = {}
with open(counted_paths_file) as f:
	for l in f:
		if l.strip() not in path_counts: path_counts[l.strip()] = 0
		path_counts[l.strip()] += 1

forbidden_paths = set()
for node in unique_nodes:
	check_side('>' + node)
	check_side('<' + node)

for path in all_paths:
	if path not in forbidden_paths:
		print(path)
	else:
		sys.stderr.write("forbade " + path + "\n")
