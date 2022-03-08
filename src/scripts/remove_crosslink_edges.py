#!/usr/bin/env python

import sys

# input graph from stdin
# output graph to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:]

def canon(left, right):
	fwstr = left + right
	bwstr = revnode(right) + revnode(left)
	if bwstr < fwstr: return (revnode(right), revnode(left))
	return (left, right)

def getone(s):
	for c in s:
		return c

def good_coverage(comparison, average):
	return comparison >= average * 0.5 and comparison <= average * 1.5

def edge_removed_oneway(kept_in_removables, fromnode, tonode):
	if fromnode not in kept_in_removables: return False
	if len(kept_in_removables[fromnode]) != 1: return False
	other = revnode(getone(kept_in_removables[fromnode]))
	if other == revnode(tonode): return False
	if other not in kept_in_removables: return False
	if len(kept_in_removables[other]) != 1: return False
	if revnode(getone(kept_in_removables[other])) != fromnode: return False
	return True

def edge_removed(kept_in_removables, fromnode, tonode):
	if edge_removed_oneway(kept_in_removables, fromnode, tonode): return True
	if edge_removed_oneway(kept_in_removables, revnode(tonode), revnode(fromnode)): return True
	return False

long_node_threshold = 100000
max_removable_coverage = 3
max_removable_len = 15000

node_lines = []
edge_lines = []
nodelens = {}
edges = {}
node_coverages = {}
edge_coverages = {}
remove_type = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		node_lines.append((parts[1], l.strip()))
		coverage = None
		for tag in parts[3:]:
			if tag[0:5] == "ll:f:":
				coverage = float(tag[5:])
		assert coverage is not None
		node_coverages[parts[1]] = coverage
		if ">" + parts[1] not in edges: edges[">" + parts[1]] = set()
		if "<" + parts[1] not in edges: edges["<" + parts[1]] = set()
		nodelens[parts[1]] = len(parts[2])
	elif parts[0] == "L":
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
		edges[fromnode].add(tonode)
		edges[revnode(tonode)].add(revnode(fromnode))
		edge_lines.append((fromnode, tonode, l.strip()))
		coverage = None
		for tag in parts[6:]:
			if tag[0:5] == "ec:i:":
				coverage = float(tag[5:])
		assert coverage is not None
		edge_coverages[canon(fromnode, tonode)] = coverage

coverage_sum = 0.0
coverage_len = 0.0

for node in nodelens:
	if nodelens[node] < long_node_threshold: continue
	coverage_len += nodelens[node]
	coverage_sum += nodelens[node] * node_coverages[node]

average_coverage = coverage_sum / coverage_len
sys.stderr.write("average coverage: " + str(average_coverage) + "\n")

assert average_coverage * 0.5 > max_removable_coverage

removed_nodes = set()

for node in node_coverages:
	if node_coverages[node] > max_removable_coverage: continue
	if nodelens[node] > max_removable_len: continue
	if len(edges[">" + node]) >= 2: continue
	if len(edges["<" + node]) >= 2: continue
	has_good_before = 1
	has_good_after = 1
	node_type = 0
	if len(edges[">" + node]) == 1:
		node_type += 1
		has_good_after = 0
		other_node = revnode(getone(edges[">" + node]))
		if good_coverage(node_coverages[other_node[1:]], average_coverage):
			for edge in edges[other_node]:
				if good_coverage(node_coverages[edge[1:]], average_coverage):
					if good_coverage(edge_coverages[canon(other_node, edge)], average_coverage):
						has_good_after += 1
	if len(edges["<" + node]) == 1:
		node_type += 1
		has_good_before = 0
		other_node = revnode(getone(edges["<" + node]))
		if good_coverage(node_coverages[other_node[1:]], average_coverage):
			for edge in edges[other_node]:
				if good_coverage(node_coverages[edge[1:]], average_coverage):
					if good_coverage(edge_coverages[canon(other_node, edge)], average_coverage):
						has_good_before += 1
	if has_good_before == 1 and has_good_after == 1:
		removed_nodes.add(node)
		remove_type[node] = node_type

kept_in_removables = {}

for fromnode in edges:
	if not good_coverage(node_coverages[fromnode[1:]], average_coverage): continue
	has_good_count = 0
	has_removable_count = 0
	for tonode in edges[fromnode]:
		if good_coverage(node_coverages[tonode[1:]], average_coverage):
			if good_coverage(edge_coverages[canon(fromnode, tonode)], average_coverage):
				has_good_count += 1
		if edge_coverages[canon(fromnode, tonode)] <= max_removable_coverage:
			has_removable_count += 1
	if has_good_count != 1 or has_good_count + has_removable_count != len(edges[fromnode]): continue
	assert fromnode not in kept_in_removables
	kept_in_removables[fromnode] = set()
	for tonode in edges[fromnode]:
		if good_coverage(node_coverages[tonode[1:]], average_coverage):
			if good_coverage(edge_coverages[canon(fromnode, tonode)], average_coverage):
				kept_in_removables[fromnode].add(tonode)

for node in removed_nodes:
	sys.stderr.write("removed node " + str(node) + " type " + str(remove_type[node]) + "\n")

for fromnode in edges:
	for tonode in edges[fromnode]:
		if edge_removed(kept_in_removables, fromnode, tonode):
			sys.stderr.write(str(fromnode) + "\t" + str(tonode) + "\n")

for line in node_lines:
	if line[0] in removed_nodes: continue
	print(line[1])

for line in edge_lines:
	if line[0][1:] in removed_nodes: continue
	if line[1][1:] in removed_nodes: continue
	if edge_removed(kept_in_removables, line[0], line[1]): continue
	print(line[2])
