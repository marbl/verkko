#!/usr/bin/env python

import sys

in_gfa = sys.argv[1]
hifi_coverage_csv = sys.argv[2]
distance_based_cuts_out = sys.argv[3]


long_node_size = 100000
min_tipnode_size = 1000
telomere_seq = "TAG"
telomere_seq_revcomp = "CTA"
telomere_end_anchorsize = 500
min_telomere_bp_count = 200
max_backtrace_length = 1000000
context_size = 10000

def iterate_deterministic(l):
	tmp = list(l)
	tmp.sort()
	for x in tmp:
		yield x

def revnode(n):
	assert len(n) >= 2
	assert n[0] == "<" or n[0] == ">"
	return (">" if n[0] == "<" else "<") + n[1:]

def find_telomeres(seq):
	if len(seq) < telomere_end_anchorsize * 2: return (False, False)
	end_count = 0
	last_match = 0
	for i in range(len(seq)-telomere_end_anchorsize, len(seq) - len(telomere_seq)):
		if seq[i:i+len(telomere_seq)] == telomere_seq:
			if i - last_match > len(telomere_seq):
				end_count += len(telomere_seq)
			else:
				end_count += i - last_match
			last_match = i
	start_count = 0
	last_match = -len(telomere_seq_revcomp)
	for i in range(0, telomere_end_anchorsize - len(telomere_seq)):
		if seq[i:i+len(telomere_seq_revcomp)] == telomere_seq_revcomp:
			if i - last_match > len(telomere_seq):
				start_count += len(telomere_seq)
			else:
				start_count += i - last_match
			last_match = i
	return (end_count > min_telomere_bp_count, start_count > min_telomere_bp_count)

def find_closest_homozygous(pos, edges, nodelens, coverage, average_coverage):
	start_coverage = coverage[pos[1:]]
	backtrace_len = 0
	visited = set()
	while backtrace_len < max_backtrace_length:
		if pos not in edges or len(edges[pos]) == 0: return None
		max_neighbor = (0, -1, 0)
		for edge in iterate_deterministic(edges[pos]):
			if edge[0] in visited: return None
			cov = 0
			if edge[0][1:] in coverage: cov = coverage[edge[0][1:]]
			if cov > max_neighbor[1]: max_neighbor = (edge[0], cov, edge[1])
		assert max_neighbor[1] >= 0
		backtrace_len += nodelens[pos[1:]] - max_neighbor[2]
		visited.add(pos)
		pos = max_neighbor[0]
		if pos[1:] in coverage and coverage[pos[1:]] >= start_coverage * 1.5 and coverage[pos[1:]] <= start_coverage * 2.5:
			return (pos, backtrace_len)
		if pos[1:] in coverage and coverage[pos[1:]] >= average_coverage * 1.5 and coverage[pos[1:]] <= average_coverage * 2.5:
			return (pos, backtrace_len)
		if pos[1:] in coverage and coverage[pos[1:]] >= average_coverage * 2.5: return None
	return None

def find_context_nodes(start, tip, backtrace_len, nodelens, edges):
	stack = [(start, 0)]
	visited = set()
	result = set()
	while len(stack) > 0:
		top = stack[-1]
		stack.pop()
		if top[0] in visited: continue
		visited.add(top[0])
		if top[0] not in edges: continue
		for edge in iterate_deterministic(edges[top[0]]):
			if edge[0] in visited: continue
			edge_start_pos = top[1] - edge[1]
			edge_end_pos = top[1] - edge[1] + nodelens[edge[0][1:]]
			if edge_end_pos >= backtrace_len - context_size and edge_start_pos <= backtrace_len + context_size:
				if edge[0] != tip and "gap" not in edge[0]:
					result.add((edge[0], edge_start_pos))
			stack.append((edge[0], edge_end_pos))
	return result

coverage = {}
with open(hifi_coverage_csv) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node": continue
		coverage[parts[0]] = float(parts[1])

not_tip = set()
maybe_tip = set()
has_telomere = set()
long_enough = set()
edges = {}
nodelens = {}

longnode_coverage = 0
longnode_sum = 0

with open(in_gfa) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			nodelens[parts[1]] = len(parts[2])
			if len(parts[2]) >= long_node_size:
				longnode_sum += len(parts[2])
				longnode_coverage += len(parts[2]) * coverage[parts[1]]
			if len(parts[2]) >= min_tipnode_size:
				maybe_tip.add(">" + parts[1])
				maybe_tip.add("<" + parts[1])
				long_enough.add(parts[1])
				(has_telomere_fw, has_telomere_bw) = find_telomeres(parts[2])
				if has_telomere_fw: has_telomere.add(">" + parts[1])
				if has_telomere_bw: has_telomere.add("<" + parts[1])
		elif parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			not_tip.add(fromnode)
			not_tip.add(revnode(tonode))
			overlap = int(parts[5][:-1])
			if fromnode not in edges: edges[fromnode] = set()
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[fromnode].add((tonode, overlap))
			edges[revnode(tonode)].add((revnode(fromnode), overlap))

average_coverage = longnode_coverage / longnode_sum

tip_nodes = set()

for tip in maybe_tip:
	if tip[1:] not in coverage: continue
	if coverage[tip[1:]] < average_coverage * 0.5 or coverage[tip[1:]] > average_coverage * 1.5: continue
	if tip[1:] not in long_enough: continue
	if tip in not_tip: continue
	if tip in has_telomere: continue
	tip_nodes.add(tip)

distance_cuts = []

for node in tip_nodes:
	hom = find_closest_homozygous(revnode(node), edges, nodelens, coverage, average_coverage)
	if not hom: continue
	context = find_context_nodes(revnode(hom[0]), node, hom[1], nodelens, edges)
	if len(context) == 0: continue
	for c in context:
		if c[1] <= hom[1] and c[1] + nodelens[c[0][1:]] >= hom[1]:
			distance_cuts.append((node, c[0], hom[1] - c[1]))

with open(distance_based_cuts_out, "w") as f:
	for cut in distance_cuts:
		f.write(cut[0] + "\t" + cut[1] + "\t" + str(cut[2]) + "\n")
