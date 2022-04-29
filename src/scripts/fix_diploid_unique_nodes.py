#!/usr/bin/env python

import sys

old_unique_nodes_file = sys.argv[1]
aln_file = sys.argv[2]
hifi_node_cov_file = sys.argv[3]
ont_node_cov_file = sys.argv[4]
graph_file = sys.argv[5]
# new uniques to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:] 

def check_side(edges, node, ont_node_covs, hifi_node_covs):
	min_ont_coverage = None
	max_ont_coverage = None
	min_hifi_coverage = None
	max_hifi_coverage = None
	for edge in edges[node]:
		if len(edges[revnode(edge)]) != 1: return False
		if min_ont_coverage is None: min_ont_coverage = ont_node_covs[edge[1:]]
		if max_ont_coverage is None: max_ont_coverage = ont_node_covs[edge[1:]]
		min_ont_coverage = min(min_ont_coverage, ont_node_covs[edge[1:]])
		max_ont_coverage = max(max_ont_coverage, ont_node_covs[edge[1:]])
		if min_hifi_coverage is None: min_hifi_coverage = hifi_node_covs[edge[1:]]
		if max_hifi_coverage is None: max_hifi_coverage = hifi_node_covs[edge[1:]]
		min_hifi_coverage = min(min_hifi_coverage, hifi_node_covs[edge[1:]])
		max_hifi_coverage = max(max_hifi_coverage, hifi_node_covs[edge[1:]])
	if max_ont_coverage > min_ont_coverage * 1.5: return False
	if max_hifi_coverage > max_hifi_coverage * 1.5: return False
	# if ont_node_covs[node[1:]] > max_ont_coverage * (len(edges[node]) + 0.5): return False
	# if ont_node_covs[node[1:]] < min_ont_coverage * (len(edges[node]) - 0.5): return False
	# if hifi_node_covs[node[1:]] > max_hifi_coverage * (len(edges[node]) + 0.5): return False
	# if hifi_node_covs[node[1:]] < min_hifi_coverage * (len(edges[node]) - 0.5): return False
	return True

def check_triplets(edges, node, triplet_coverages):
	copycount = len(edges[node])
	if node not in triplet_coverages: return False
	if len(triplet_coverages[node]) < copycount: return False
	counts = [triplet_coverages[node][key] for key in triplet_coverages[node]]
	counts.sort(key = lambda x: -x)
	if counts[0] > counts[copycount-1] * 2: return False
	out_covered = set()
	in_covered = set()
	triplets = list(triplet_coverages[node])
	triplets.sort(key = lambda x: -triplet_coverages[node][x])
	for i in range(0, copycount):
		out_covered.add(triplets[i][2])
		in_covered.add(triplets[i][0])
	if len(in_covered) != len(edges[revnode(node)]): return False
	if len(out_covered) != len(edges[node]): return False
	if len(triplet_coverages[node]) > copycount and counts[copycount-1] < counts[copycount] * 2: return False
	return True

unique_nodes = set()
with open(old_unique_nodes_file) as f:
	for l in f:
		unique_nodes.add(l.strip())

ont_node_covs = {}
with open(ont_node_cov_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'node': continue
		assert parts[0] not in ont_node_covs
		ont_node_covs[parts[0]] = float(parts[2])

hifi_node_covs = {}
with open(hifi_node_cov_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'node': continue
		assert parts[0] not in hifi_node_covs
		hifi_node_covs[parts[0]] = float(parts[1])

nodes = set()
edges = {}
check_node_queue = []
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if parts[1] not in ont_node_covs: ont_node_covs[parts[1]] = 0
			check_node_queue.append(parts[1])
		elif parts[0] == 'L':
			fromnode = ('>' if parts[2] == '+' else '<') + parts[1]
			tonode = ('>' if parts[4] == '+' else '<') + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))

triplet_coverages = {}
with open(aln_file) as f:
	for l in f:
		path = l.strip().split('\t')[5].replace('>', '\t>').replace('<', '\t<').strip().split('\t')
		for i in range(1, len(path)-1):
			triplet = (path[i-1], path[i], path[i+1])
			if path[i][0] == "<":
				triplet = (revnode(path[i+1]), revnode(path[i]), revnode(path[i-1]))
			assert triplet[1][0] == ">"
			if triplet[1] not in triplet_coverages: triplet_coverages[triplet[1]] = {}
			if triplet not in triplet_coverages[triplet[1]]: triplet_coverages[triplet[1]][triplet] = 0
			triplet_coverages[triplet[1]][triplet] += 1

while len(check_node_queue) > 0:
	node = check_node_queue.pop()
	if node in unique_nodes: continue
	if '>' + node not in edges: continue
	if '<' + node not in edges: continue
	if len(edges['>' + node]) != len(edges['<' + node]): continue
	if len(edges['>' + node]) < 2: continue
	fw_has_unique = False
	bw_has_unique = False
	fw_all_unique = True
	bw_all_unique = True
	for edge in edges['>' + node]:
		if edge[1:] in unique_nodes:
			fw_has_unique = True
		else:
			fw_all_unique = False
	for edge in edges['<' + node]:
		if edge[1:] in unique_nodes:
			bw_has_unique = True
		else:
			bw_all_unique = False
	if fw_all_unique and bw_all_unique: continue
	if not fw_has_unique and not bw_has_unique: continue
	assert fw_has_unique or bw_has_unique
	assert not fw_all_unique or not bw_all_unique
	fw_valid = check_side(edges, ">" + node, ont_node_covs, hifi_node_covs)
	if not fw_valid: continue
	bw_valid = check_side(edges, "<" + node, ont_node_covs, hifi_node_covs)
	if not bw_valid: continue
	triplets_valid = check_triplets(edges, ">" + node, triplet_coverages)
	if not triplets_valid: continue
	for edge in edges[">" + node]:
		unique_nodes.add(edge[1:])
		if edge in edges:
			for edge2 in edges[edge]:
				check_node_queue.append(edge2[1:])
	for edge in edges["<" + node]:
		unique_nodes.add(edge[1:])
		if edge in edges:
			for edge2 in edges[edge]:
				check_node_queue.append(edge2[1:])

for node in unique_nodes: print(node)

