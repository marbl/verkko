#!/usr/bin/env python

import sys
import graph_functions as gf

old_unique_nodes_file = sys.argv[1]
coverage_file = sys.argv[2]
# graph from stdin
# new uniques to stdout

def get_reachable_uniques(start, edges, uniques):
	checked = set()
	stack = [node]
	result = set()
	while len(stack) > 0:
		top = stack[-1]
		stack.pop()
		if top in checked: continue
		checked.add(top)
		if top in uniques:
			result.add(top)
			continue
		if ">" + top in edges:
			for edge in edges[">" + top]:
				stack.append(edge[1:])
		if "<" + top in edges:
			for edge in edges["<" + top]:
				stack.append(edge[1:])
	return result

uniques = set()
with open(old_unique_nodes_file) as f:
	for l in f:
		uniques.add(l.strip())

coverages = {}
with open(coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node":
			assert(parts[1] == "coverage" and parts[2] == "length")
			continue
		coverages[parts[0]] = float(parts[1])

nodelens = {}
edges = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] == "S":
		nodelens[parts[1]] = len(parts[2])
	if parts[0] == "L":
		fromnode = (">" if parts[2] == "+" else "<") + parts[1]
		tonode = (">" if parts[4] == "+" else "<") + parts[3]
		if fromnode not in edges: edges[fromnode] = set()
		edges[fromnode].add(tonode)
		if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
		edges[gf.revnode(tonode)].add(gf.revnode(fromnode))

coverage_sum = 0.0
coverage_len = 0

for node in uniques:
	assert node in nodelens
	if node not in coverages: continue
	coverage_sum += coverages[node] * nodelens[node]
	coverage_len += nodelens[node]

avg_coverage = coverage_sum / coverage_len if float(coverage_len) > 0.0 else 0.0

new_uniques = set()

for node in coverages:
	if node in uniques: continue
	if coverages[node] < avg_coverage * .5 or coverages[node] > avg_coverage * 1.5: continue
	if ">" + node not in edges: continue
	if "<" + node not in edges: continue
	if len(edges[">" + node]) != 1: continue
	if len(edges["<" + node]) != 1: continue
	if gf.getone(edges[">" + node]) != gf.getone(edges["<" + node]): continue
#	if len(get_reachable_uniques(node, edges, uniques)) != 2: continue
	new_uniques.add(node)

for node in new_uniques:
	assert node not in uniques
	uniques.add(node)
	sys.stderr.write(node + "\n")

for node in uniques:
	print(node)
