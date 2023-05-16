#!/usr/bin/env python

import sys

input_resolved_gfa = sys.argv[1]
input_nodemap = sys.argv[2]
gap_prefix = sys.argv[3]
# output to stdout

max_acyclic_continuation_size = 200000

def get_continuation_acyclic_size(start_node, edges, node_lens):
	visited = set()
	check_stack = []
	if start_node not in edges: return 0
	for edge in edges[start_node]:
		check_stack.append(edge[1:])
	result = 0
	visited.add(start_node[1:])
	while len(check_stack) > 0:
		top = check_stack[-1]
		check_stack.pop()
		if top in visited: continue
		visited.add(top)
		result += node_lens[top]
		if ">" + top in edges:
			for edge in edges[">" + top]:
				check_stack.append(edge[1:])
		if "<" + top in edges:
			for edge in edges["<" + top]:
				check_stack.append(edge[1:])
	return result

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(n1, n2):
	fw = n1 + n2
	bw = revnode(n2) + revnode(n1)
	if bw < fw: return (revnode(n2), revnode(n1))
	return (n1, n2)

def getone(s):
	for c in s:
		return c

def reduce_to_leaves(path, nodemap):
	result = list(path)
	while True:
		changed = False
		new_result = []
		for node in result:
			if node[1:] not in nodemap:
				new_result.append(node)
			else:
				add = nodemap[node[1:]]
				if node[0] == "<":
					add = [revnode(n) for n in add[::-1]]
				new_result += add
				changed = True
		if not changed: break
		result = new_result
	return result

def find(parent, key):
	while parent[key] != parent[parent[key]]: parent[key] = parent[parent[key]]
	return parent[key]

def merge(parent, left, right):
	left = find(parent, left)
	right = find(parent, right)
	assert parent[left] == left
	assert parent[right] == right
	parent[right] = left

nodemap = {}
with open(input_nodemap) as f:
	for l in f:
		parts = l.strip().split("\t")
		keynode = parts[0]
		splitnodes = parts[1].split(":")[0].replace(">", "\t>").replace("<", "\t<").strip().split("\t")
		assert keynode not in nodemap or nodemap[keynode] == splitnodes
		nodemap[keynode] = splitnodes

node_parent = {}
final_graph_nodes = set()
node_lens = {}
edges = {}
with open(input_resolved_gfa) as f:
	for l in f:
		print(l.strip())
		parts = l.strip().split("\t")
		if parts[0] == "S":
			if parts[1] not in node_parent: node_parent[parts[1]] = parts[1]
			node_lens[parts[1]] = len(parts[2])
			final_graph_nodes.add(parts[1])
		if parts[0] == "L":
			if parts[1] not in node_parent: node_parent[parts[1]] = parts[1]
			if parts[3] not in node_parent: node_parent[parts[3]] = parts[1]
			merge(node_parent, parts[1], parts[3])
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))

leaf_belongs_to_two_nodes = {}
for node in final_graph_nodes:
	leaves = reduce_to_leaves([">" + node], nodemap)
	for leaf in leaves:
		if leaf[1:] not in leaf_belongs_to_two_nodes:
			leaf_belongs_to_two_nodes[leaf[1:]] = (leaf[0] + node, None)
		elif leaf[1:] in leaf_belongs_to_two_nodes and type(leaf_belongs_to_two_nodes[leaf[1:]]) is tuple and leaf_belongs_to_two_nodes[leaf[1:]][1] is None and leaf_belongs_to_two_nodes[leaf[1:]][0][1:] != node:
			leaf_belongs_to_two_nodes[leaf[1:]] = (leaf_belongs_to_two_nodes[leaf[1:]][0], leaf[0] + node)
		else:
			leaf_belongs_to_two_nodes[leaf[1:]] = None

possible_edges = {}
for node in final_graph_nodes:
	leaves = reduce_to_leaves([">" + node], nodemap)
	last_match = None
	for leaf in leaves:
		assert leaf[1:] in leaf_belongs_to_two_nodes
		if leaf_belongs_to_two_nodes[leaf[1:]] is None: continue
		assert type(leaf_belongs_to_two_nodes[leaf[1:]]) is tuple
		assert len(leaf_belongs_to_two_nodes[leaf[1:]]) == 2
		if leaf_belongs_to_two_nodes[leaf[1:]][1] is None: continue
		assert leaf_belongs_to_two_nodes[leaf[1:]][0][1:] == node or leaf_belongs_to_two_nodes[leaf[1:]][1][1:] == node
		assert leaf_belongs_to_two_nodes[leaf[1:]][0][1:] != node or leaf_belongs_to_two_nodes[leaf[1:]][1][1:] != node
		if leaf_belongs_to_two_nodes[leaf[1:]][0][1:] == node:
			this_orient = leaf_belongs_to_two_nodes[leaf[1:]][0][0]
			other_node = leaf_belongs_to_two_nodes[leaf[1:]][1]
		else:
			this_orient = leaf_belongs_to_two_nodes[leaf[1:]][1][0]
			other_node = leaf_belongs_to_two_nodes[leaf[1:]][0]
		if this_orient == "<":
			other_node = revnode(other_node)
		if last_match is None:
			last_match = other_node
			continue
		if other_node == last_match: continue
		if find(node_parent, last_match[1:]) == find(node_parent, node) or find(node_parent, other_node[1:]) == find(node_parent, node):
			if find(node_parent, last_match[1:]) != find(node_parent, other_node[1:]):
				if get_continuation_acyclic_size(last_match, edges, node_lens) < max_acyclic_continuation_size:
					if get_continuation_acyclic_size(revnode(other_node), edges, node_lens) < max_acyclic_continuation_size:
						if last_match not in possible_edges: possible_edges[last_match] = set()
						possible_edges[last_match].add(other_node)
						if revnode(other_node) not in possible_edges: possible_edges[revnode(other_node)] = set()
						possible_edges[revnode(other_node)].add(revnode(last_match))
						sys.stderr.write("potential edge " + last_match + " to " + other_node + "\n")
		last_match = other_node

next_gap_num = 0
added_gaps = set()
gap_sequence = "N"*20000
for edge in possible_edges:
	if len(possible_edges[edge]) != 1: continue
	other_side = getone(possible_edges[edge])
	if len(possible_edges[revnode(other_side)]) != 1: continue
	assert getone(possible_edges[revnode(other_side)]) == revnode(edge)
	key = canon(edge, other_side)
	if key in added_gaps: continue
	added_gaps.add(key)
	gap_name = gap_prefix + str(next_gap_num)
	next_gap_num += 1
	sys.stderr.write("add edge from " + edge + " to " + other_side + "\n")
	print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + other_side[1:] + "\t" + ("+" if other_side[0] == ">" else "-") + "\t0M")
