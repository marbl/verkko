#!/usr/bin/python

import sys

old_unique_nodes_file = sys.argv[1]
node_cov_file = sys.argv[2]
graph_file = sys.argv[3]
# new uniques to stdout

def revnode(n):
	assert len(n) >= 2
	assert n[0] == '>' or n[0] == '<'
	return ('>' if n[0] == '<' else '<') + n[1:] 

unique_nodes = set()
with open(old_unique_nodes_file) as f:
	for l in f:
		unique_nodes.add(l.strip())

node_covs = {}
with open(node_cov_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'node': continue
		assert parts[0] not in node_covs
		node_covs[parts[0]] = float(parts[2])

nodes = set()
edges = {}
check_node_queue = []
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if parts[1] not in node_covs: node_covs[parts[1]] = 0
			check_node_queue.append(parts[1])
		elif parts[0] == 'L':
			fromnode = ('>' if parts[2] == '+' else '<') + parts[1]
			tonode = ('>' if parts[4] == '+' else '<') + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			edges[fromnode].add(tonode)
			if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
			edges[revnode(tonode)].add(revnode(fromnode))

while len(check_node_queue) > 0:
	node = check_node_queue.pop()
	if node in unique_nodes: continue
	if '>' + node not in edges: continue
	if '<' + node not in edges: continue
	if len(edges['>' + node]) != len(edges['<' + node]): continue
	if len(edges['>' + node]) < 2: continue
	fw_all_unique = True
	bw_all_unique = True
	for edge in edges['>' + node]:
		if edge[1:] not in unique_nodes:
			fw_all_unique = False
			break
	for edge in edges['<' + node]:
		if edge[1:] not in unique_nodes:
			bw_all_unique = False
	if fw_all_unique and bw_all_unique: continue
	if not fw_all_unique and not bw_all_unique: continue
	assert fw_all_unique or bw_all_unique
	assert not fw_all_unique or not bw_all_unique
	if fw_all_unique:
		assert fw_all_unique
		assert not bw_all_unique
		unique_side = '>' + node
		not_unique_side = '<' + node
	else:
		assert not fw_all_unique
		assert bw_all_unique
		unique_side = '<' + node
		not_unique_side = '>' + node
	unique_side_valid = True
	for edge in edges[unique_side]:
		if len(edges[revnode(edge)]) != 1:
			unique_side_valid = False
			break
	if not unique_side_valid: continue
	not_unique_side_valid = True
	min_coverage = None
	max_coverage = None
	for edge in edges[not_unique_side]:
		if len(edges[revnode(edge)]) != 1:
			not_unique_side_valid = False
			break
		if min_coverage is None: min_coverage = node_covs[edge[1:]]
		if max_coverage is None: max_coverage = node_covs[edge[1:]]
		min_coverage = min(min_coverage, node_covs[edge[1:]])
		max_coverage = max(max_coverage, node_covs[edge[1:]])
	if not not_unique_side_valid: continue
	if max_coverage > min_coverage * 1.5: continue
	for edge in edges[not_unique_side]:
		if edge[1:] in unique_nodes: continue
		unique_nodes.add(edge[1:])
		if edge in edges:
			for edge2 in edges[edge]:
				check_node_queue.append(edge2[1:])

for node in unique_nodes: print(node)

