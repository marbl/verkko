#!/usr/bin/env python

import sys
import graph_functions as gf

graph_file = sys.argv[1]
longest_n = int(sys.argv[2])
# paths from stdin
# graph to stdout

max_unroll_length = 200000
cov_file = sys.argv[3]

def iscanon(left, right):
        fwstr = left + right
        bwstr = right + left
        if bwstr < fwstr: return False
        return True

coverage = {}
min_len = max_unroll_length / 2
long_coverage_len_sum = 0
long_coverage_cov_sum = 0
coverage_len_sum = 0
coverage_cov_sum = 0
with open(cov_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if "node" in parts[0]: continue
		coverage[parts[0]] = float(parts[1])
		nodelen = int(parts[2])
		coverage_len_sum += nodelen
		coverage_cov_sum += nodelen * coverage[parts[0]]
		if nodelen > min_len:
			long_coverage_len_sum += nodelen
			long_coverage_cov_sum += nodelen * coverage[parts[0]]
avg_coverage = float(coverage_cov_sum) / float(coverage_len_sum)
if long_coverage_len_sum > 0:
	avg_coverage = float(long_coverage_cov_sum) / float(long_coverage_len_sum)

nodeseqs = {}
edges = {}
overlaps = {}
max_overlap = {}
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
			if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
			edges[gf.revnode(tonode)].add(gf.revnode(fromnode))
			overlaps[(fromnode, tonode)] = parts[5]
			overlaps[(gf.revnode(tonode), gf.revnode(fromnode))] = parts[5]

			if parts[1] != parts[3] and iscanon(fromnode, tonode):
				if fromnode not in max_overlap: max_overlap[fromnode] = int(parts[5][:-1])
				max_overlap[fromnode] = max(max_overlap[fromnode], int(parts[5][:-1]))
				if tonode not in max_overlap: max_overlap[tonode] = int(parts[5][:-1])
				max_overlap[tonode] = max(max_overlap[tonode], int(parts[5][:-1]))

tip_loops = set()
counts_per_tiploop = {}

for node in nodeseqs:
	if ">" + node not in edges: continue
	if "<" + node not in edges: continue

	if len(edges[">" + node]) > 2 or len(edges["<" + node]) > 2: continue
	if ">" + node not in edges[">" + node]: continue
	if len(edges[">" + node]) == 2 and len(edges["<" + node]) == 2: continue
	if len(nodeseqs[node]) > max_unroll_length: continue
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
	if node in coverage and coverage[node] < 1.1*avg_coverage: continue  # don't unroll nodes that seem to have insufficient coverage to be a loop
	assert count >= 1
	assert len(edges["<" + node]) >= 1 and len(edges["<" + node]) <= 2
	assert len(edges[">" + node]) >= 1 and len(edges[">" + node]) <= 2
	if len(edges["<" + node]) == 2:
		for edge in edges["<" + node]:
			if edge[1:] == node: continue
			assert ">" + node in edges[gf.revnode(edge)]
			edges[gf.revnode(edge)].remove(">" + node)
			edges[gf.revnode(edge)].add(">unroll_" + node + "_1")
			edges["<unroll_" + node + "_1"] = set()
			edges["<unroll_" + node + "_1"].add(edge)
			overlaps[(gf.revnode(edge), ">unroll_" + node + "_1")] = overlaps[(gf.revnode(edge), ">" + node)]
			overlaps[("<unroll_" + node + "_1", edge)] = overlaps[(gf.revnode(edge), ">" + node)]
	self_overlap = overlaps[(">" + node, ">" + node)]
	for i in range(0, count):
		nodeseqs["unroll_" + node + "_" + str(i+1)] = nodeseqs[node]
		edges[">" + "unroll_" + node + "_" + str(i+1)] = set()
		if i > 0: edges["<" + "unroll_" + node + "_" + str(i+1)] = set()
	for i in range(1, count):
		edges[">" + "unroll_" + node + "_" + str(i)].add(">" +  "unroll_" + node + "_" + str(i+1))
		overlaps[(">" + "unroll_" + node + "_" + str(i), ">" +  "unroll_" + node + "_" + str(i+1))] = self_overlap
	if len(edges[">" + node]) == 2:
		for edge in edges[">" + node]:
			if edge[1:] == node: continue
			assert "<" + node in edges[gf.revnode(edge)]
			edges[gf.revnode(edge)].remove("<" + node)
			edges[gf.revnode(edge)].add("<unroll_" + node + "_" + str(count))
			assert ">unroll_" + node + "_" + str(count) in edges
			edges[">unroll_" + node + "_"  + str(count)].add(edge)
			overlaps[(gf.revnode(edge), "<unroll_" + node + "_" + str(count))] = overlaps[(gf.revnode(edge), "<" + node)]
			overlaps[(">unroll_" + node + "_" + str(count), edge)] = overlaps[(gf.revnode(edge), "<" + node)]
	del edges[">" + node]
	del edges["<" + node]
	del nodeseqs[node]
	unrolled_count += 1
	sys.stderr.write("unroll " + str(node) + "\n")

sys.stderr.write("unrolled " + str(unrolled_count) + " looptips\n")

for node in tip_loops:
	if node not in nodeseqs: continue

	nodelen=len(nodeseqs[node])
	if ">"+node in max_overlap: nodelen -= max_overlap[">"+node]
	if "<"+node in max_overlap: nodelen -= max_overlap["<"+node]
	if nodelen <= 0: nodelen = 1

	# if a node is short enough (we tolerate a bit longer than the unrolling length and doesn't look like it's high coverage enough to be unrolled then we drop it's loop edge and hope we can reconnect
	if node in coverage and coverage[node] < 1.1*avg_coverage and nodelen <= 1.5*max_unroll_length:
		sys.stderr.write("%s removed edges\n"%(node))
		del edges[">" + node]
		del edges["<" + node]

nodenames = list(nodeseqs)
nodenames.sort()

for node in nodenames:
	print("S\t" + node + "\t" + nodeseqs[node])

edgenames = list(edges)
edgenames.sort()

for edge in edgenames:
	targets = list(edges[edge])
	targets.sort()
	for target in targets:
		print("L\t" + edge[1:] + "\t" + ("+" if edge[0] == ">" else "-") + "\t" + target[1:] + "\t" + ("+" if target[0] == ">" else "-") + "\t" + overlaps[(edge, target)])
