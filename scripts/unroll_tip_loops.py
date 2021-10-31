#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
longest_n = int(sys.argv[2])
# paths from stdin
# graph to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

nodeseqs = {}
edges = {}
overlaps = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodeseqs[parts[1]] = parts[2]
		elif parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))
			overlaps[(fromnode, tonode)] = parts[5]
			overlaps[(revnode(tonode), revnode(fromnode))] = parts[5]

tip_loops = set()
counts_per_tiploop = {}

for node in nodeseqs:
	if ">" + node not in edges: continue
	if "<" + node not in edges: continue
	if len(edges[">" + node]) > 2 or len(edges["<" + node]) > 2: continue
	if ">" + node not in edges[">" + node]: continue
	if len(edges[">" + node]) == 2 and len(edges["<" + node]) == 2: continue
	tip_loops.add(node)
	counts_per_tiploop[node] = []

for l in sys.stdin:
	path = l.replace('>', '\t>').replace('<', '\t<').strip().split('\t')
	counts = {}
	for node in path:
		if node[1:] not in tip_loops: continue
		if node[1:] not in counts: counts[node[1:]] = 0
		counts[node[1:]] += 1
	for node in counts:
		assert node in counts_per_tiploop
		counts_per_tiploop[node].append(counts[node])

unrolled_count = 0

for node in counts_per_tiploop:
	if len(counts_per_tiploop[node]) < longest_n: continue
	counts_per_tiploop[node].sort()
	assert len(counts_per_tiploop[node]) >= longest_n
	count = counts_per_tiploop[node][-longest_n]
	assert count >= 1
	assert len(edges["<" + node]) >= 1 and len(edges["<" + node]) <= 2
	assert len(edges[">" + node]) >= 1 and len(edges[">" + node]) <= 2
	if len(edges["<" + node]) == 2:
		for edge in edges["<" + node]:
			if edge[1:] == node: continue
			assert ">" + node in edges[revnode(edge)]
			edges[revnode(edge)].remove(">" + node)
			edges[revnode(edge)].add(">unroll_" + node + "_1")
			overlaps[(revnode(edge), ">unroll_" + node + "_1")] = overlaps[(revnode(edge), ">" + node)]
			overlaps[("<unroll_" + node + "_1", edge)] = overlaps[(revnode(edge), ">" + node)]
	self_overlap = overlaps[(">" + node, ">" + node)]
	for i in range(0, count):
		nodeseqs["unroll_" + node + "_" + str(i+1)] = nodeseqs[node]
		edges[">" + "unroll_" + node + "_" + str(i+1)] = set()
		edges["<" + "unroll_" + node + "_" + str(i+1)] = set()
	for i in range(1, count):
		edges[">" + "unroll_" + node + "_" + str(i)].add(">" +  "unroll_" + node + "_" + str(i+1))
		overlaps[(">" + "unroll_" + node + "_" + str(i), ">" +  "unroll_" + node + "_" + str(i+1))] = self_overlap
	if len(edges[">" + node]) == 2:
		for edge in edges[">" + node]:
			if edge[1:] == node: continue
			#assert "<" + node in edges[revnode(edge)]
			if "<" + node not in edges[revnode(edge)]:
				print("WARN: prevented assertion check for <", node)
				continue
			edges[revnode(edge)].remove("<" + node)
			edges[revnode(edge)].add("<unroll_" + node + "_" + str(count))
			overlaps[(revnode(edge), "<unroll_" + node + "_" + str(count))] = overlaps[(revnode(edge), "<" + node)]
			overlaps[(">unroll_" + node + "_" + str(count), edge)] = overlaps[(revnode(edge), "<" + node)]
	del edges[">" + node]
	del edges["<" + node]
	del nodeseqs[node]
	unrolled_count += 1
	sys.stderr.write("unroll " + str(node) + "\n")

sys.stderr.write("unrolled " + str(unrolled_count) + " looptips\n")

for node in nodeseqs:
	print("S\t" + node + "\t" + nodeseqs[node])
for edge in edges:
	for target in edges[edge]:
		print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + target[1:] + "\t" + ("+" if target[0] == ">" else "-") + "\t" + overlaps[(edge, target)])
