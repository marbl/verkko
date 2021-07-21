#!/usr/bin/python

import sys

graph_file = sys.argv[1]
alignment_file = sys.argv[2]
long_node_threshold = int(sys.argv[3])
solid_edge_threshold = int(sys.argv[4])
# unique nodes to stdout

long_node_neighborhood_size = long_node_threshold

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(left, right):
	fwstr = left + right
	bwstr = revnode(right) + revnode(left)
	if bwstr < fwstr: return (revnode(right), revnode(left))
	return (left, right)

def canontip(left, right):
	fwstr = left + right
	bwstr = right + left
	if bwstr < fwstr: return (right, left)
	return (left, right)

def find(parent, key):
	assert key in parent
	while parent[key] != parent[parent[key]]:
		parent[key] = parent[parent[key]]
	return parent[key]

def merge(parent, rank, left, right):
	left = find(parent, left)
	right = find(parent, right)
	assert left in rank
	assert right in rank
	if rank[left] < rank[right]: (left, right) = (right, left)
	parent[right] = left
	if rank[left] == rank[right]: rank[right] += 1

def getone(s):
	assert len(s) == 1
	for i in s:
		return i

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013
# fig. 5
def find_bubble_end(edges, s):
	if s not in edges: return None
	if len(edges[s]) < 2: return None
	S = [s]
	visited = set()
	seen = set()
	seen.add(s)
	while len(S) > 0:
		v = S.pop()
		assert v in seen
		seen.remove(v)
		assert v not in visited
		visited.add(v)
		if v not in edges: return None
		if len(edges[v]) == 0: return None
		for u in edges[v]:
			if u[1:] == v[1:]: return None
			if revnode(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert revnode(u) in edges
			assert len(edges[revnode(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in edges[revnode(u)]:
				parent = revnode(parent_edge)
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == getone(seen):
			t = S.pop()
			if t in edges:
				for edge in edges[t]:
					if edge == s: return None
			return t
	return None

parent = {}
rank = {}

existing_nodes = set()
long_nodes = set()
canon_edges = set()
edges = {}
nodelens = {}
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			existing_nodes.add(parts[1])
			nodelens[parts[1]] = len(parts[2])
			if len(parts[2]) > long_node_threshold:
				long_nodes.add(parts[1])
		elif parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromnode not in edges: edges[fromnode] = set()
			if tonode not in edges: edges[tonode] = set()
			edges[fromnode].add(revnode(tonode))
			edges[tonode].add(revnode(fromnode))
			canon_edges.add(canontip(fromnode, tonode))

for node in nodelens:
	parent[">" + node] = ">" + node
	rank[">" + node] = 1
	parent["<" + node] = "<" + node
	rank["<" + node] = 1
	if node not in long_nodes: merge(parent, rank, ">" + node, "<" + node)

for edge in canon_edges:
	merge(parent, rank, edge[0], edge[1])

cluster_edge_nodes = {}
longnode_coverage = {}

for node in long_nodes:
	longnode_coverage[">" + node] = 0.0
	longnode_coverage["<" + node] = 0.0
	key = find(parent, ">" + node)
	if key not in cluster_edge_nodes: cluster_edge_nodes[key] = set()
	cluster_edge_nodes[key].add(">" + node)
	key = find(parent, "<" + node)
	if key not in cluster_edge_nodes: cluster_edge_nodes[key] = set()
	cluster_edge_nodes[key].add("<" + node)

node_coverage = {}
for node in nodelens:
	node_coverage[node] = 0.0

out_paths = {}

with open(alignment_file) as f:
	for l in f:
		parts = l.split('\t')
		left_clip = int(parts[7])
		right_clip = int(parts[6]) - int(parts[8])
		path = parts[5].replace(">", "\t>").replace("<", "\t<").strip().split('\t')
		part_paths = []
		last_break = 0
		for i in range(0, len(path)):
			if path[i][1:] not in existing_nodes:
				if i > last_break: part_paths.append(path[last_break:i])
				last_break = i+1
		if last_break < len(path): part_paths.append(path[last_break:])
		for part_path in part_paths:
			assert len(part_path) >= 1
			for i in range(0, len(part_path)):
				if part_path[i] not in out_paths: out_paths[part_path[i]] = set()
				if revnode(part_path[i]) not in out_paths: out_paths[revnode(part_path[i])] = set()
				if i < len(part_path)-1: out_paths[part_path[i]].add(tuple(part_path[i+1:]))
				if i > 0: out_paths[revnode(part_path[i])].add(tuple(revnode(n) for n in part_path[0:i][::-1]))
		assert len(path) >= 1
		while len(path) > 0 and path[0][1:] not in nodelens:
			path = path[1:]
			left_clip = 0
		while len(path) > 0 and path[-1][1:] not in nodelens:
			path = path[:-1]
			right_clip = 0
		if len(path) == 0: continue
		if len(path) == 1:
			# should be true, except for a graphaligner bug in self-loop nodes causing strange alignment paths
			# so comment for now
			# assert left_clip + right_clip < nodelens[path[0][1:]]
			if left_clip + right_clip >= nodelens[path[0][1:]]: continue
			node_coverage[path[0][1:]] += float(nodelens[path[0][1:]] - left_clip - right_clip) / float(nodelens[path[0][1:]])
			if path[0][1:] in long_nodes:
				node_start = left_clip
				node_end = nodelens[path[0][1:]] - right_clip
				assert node_end > node_start
				if node_start < long_node_neighborhood_size: longnode_coverage[revnode(path[0])] += float(min(long_node_neighborhood_size, node_end) - node_start) / float(long_node_neighborhood_size)
				if node_end > nodelens[path[0][1:]] - long_node_neighborhood_size: longnode_coverage[path[0]] += float(node_end - max(nodelens[path[0][1:]] - long_node_neighborhood_size, node_start)) / float(long_node_neighborhood_size)
			continue
		assert left_clip < nodelens[path[0][1:]]
		assert right_clip < nodelens[path[-1][1:]]
		assert path[0][1:] in existing_nodes
		assert path[-1][1:] in existing_nodes
		node_coverage[path[0][1:]] += float(nodelens[path[0][1:]] - left_clip) / float(nodelens[path[0][1:]])
		node_coverage[path[-1][1:]] += float(nodelens[path[-1][1:]] - right_clip) / float(nodelens[path[-1][1:]])
		for node in path[1:-1]:
			if node[1:] not in existing_nodes: continue
			node_coverage[node[1:]] += 1
			if node[1:] in long_nodes:
				longnode_coverage[">" + node[1:]] += 1
				longnode_coverage["<" + node[1:]] += 1
		if path[0][1:] in long_nodes:
			node_start = left_clip
			node_end = nodelens[path[0][1:]]
			assert node_end > node_start
			if node_start < long_node_neighborhood_size: longnode_coverage[revnode(path[0])] += float(min(long_node_neighborhood_size, node_end) - node_start) / float(long_node_neighborhood_size)
			if node_end > nodelens[path[0][1:]] - long_node_neighborhood_size: longnode_coverage[path[0]] += float(node_end - max(nodelens[path[0][1:]] - long_node_neighborhood_size, node_start)) / float(long_node_neighborhood_size)
		if path[-1][1:] in long_nodes:
			node_start = 0
			node_end = nodelens[path[-1][1:]] - right_clip
			assert node_end > node_start
			if node_start < long_node_neighborhood_size: longnode_coverage[revnode(path[-1])] += float(min(long_node_neighborhood_size, node_end) - node_start) / float(long_node_neighborhood_size)
			if node_end > nodelens[path[-1][1:]] - long_node_neighborhood_size: longnode_coverage[path[-1]] += float(node_end - max(nodelens[path[-1][1:]] - long_node_neighborhood_size, node_start)) / float(long_node_neighborhood_size)

coverage_sum = 0.0
coverage_count = 0.0
for node in long_nodes:
	coverage_sum += node_coverage[node] * nodelens[node]
	coverage_count += nodelens[node]
global_average_coverage = float(coverage_sum) / float(coverage_count)

normalized_node_coverage = {}
for node in nodelens:
	if node in long_nodes:
		normalized_node_coverage[node] = 1.0
		continue
	key = find(parent, ">" + node)
	compare_nodes = set()
	if key in cluster_edge_nodes: compare_nodes = cluster_edge_nodes[key]
	if len(compare_nodes) == 0:
		compare_coverage = global_average_coverage
	else:
		compare_coverage = 0.0
		for end in compare_nodes:
			compare_coverage += float(longnode_coverage[end])
		compare_coverage /= float(len(compare_nodes))
	normalized_node_coverage[node] = node_coverage[node] / compare_coverage

for n in out_paths:
	if n[0] != ">": continue
	if revnode(n) not in out_paths: continue
	if n[1:] not in normalized_node_coverage: continue
	if normalized_node_coverage[n[1:]] < 0.5 and node_coverage[n[1:]] / global_average_coverage < 0.5: continue
	if normalized_node_coverage[n[1:]] > 1.5 and node_coverage[n[1:]] / global_average_coverage > 1.5: continue
	inconsistent = False
	fw_paths = list(out_paths[n])
	fw_paths.sort(key=lambda x: len(x))
	for path in fw_paths:
		if path != fw_paths[-1][0:len(path)]: inconsistent = True
	bw_paths = list(out_paths[revnode(n)])
	bw_paths.sort(key=lambda x: len(x))
	for path in bw_paths:
		if path != bw_paths[-1][0:len(path)]: inconsistent = True
	if inconsistent: continue
	print(n[1:])

del out_paths

chain_parent = {}
chain_rank = {}
for node in nodelens:
	chain_parent[node] = node
	chain_rank[node] = 0

for node in nodelens:
	bubble_end = find_bubble_end(edges, ">" + node)
	if bubble_end:
		merge(chain_parent, chain_rank, node, bubble_end[1:])
	bubble_end = find_bubble_end(edges, "<" + node)
	if bubble_end:
		merge(chain_parent, chain_rank, node, bubble_end[1:])

chain_coverage_sum = {}
chain_coverage_count = {}

for node in nodelens:
	chain = find(chain_parent, node)
	if chain not in chain_coverage_count:
		assert chain not in chain_coverage_sum
		chain_coverage_count[chain] = 0.0
		chain_coverage_sum[chain] = 0.0
	if nodelens[node] < 20000: continue
	chain_coverage_count[chain] += nodelens[node]
	chain_coverage_sum[chain] += nodelens[node] * normalized_node_coverage[node]

unique_chains = set()
for node in nodelens:
	chain = find(chain_parent, node)
	assert chain in chain_coverage_count
	assert chain in chain_coverage_sum
	if chain_coverage_count[chain] == 0: continue
	if chain_coverage_sum[chain] == 0: continue
	chain_coverage = float(chain_coverage_sum[chain]) / float(chain_coverage_count[chain])
	chain_length = float(chain_coverage_sum[chain])
	unique_chain = False
	if chain_length > 35000 and chain_coverage > 0.9 and chain_coverage < 1.1: unique_chain = True
	if chain_length > 50000 and chain_coverage > 0.8 and chain_coverage < 1.2: unique_chain = True
	if chain_length > 75000 and chain_coverage > 0.75 and chain_coverage < 1.25: unique_chain = True
	if chain_length > 100000 and chain_coverage > 0.7 and chain_coverage < 1.3: unique_chain = True
	if chain_length > 200000 and chain_coverage > 0.6 and chain_coverage < 1.4: unique_chain = True
	if unique_chain: unique_chains.add(chain)

for node in nodelens:
	chain = find(chain_parent, node)
	if chain not in unique_chains: continue
	if find_bubble_end(edges, ">" + node):
		for edge in edges[">" + node]:
			if find(chain_parent, edge[1:]) in unique_chains:
				unique_chains.remove(chain)
				break
	if chain not in unique_chains: continue
	if find_bubble_end(edges, "<" + node):
		for edge in edges["<" + node]:
			if find(chain_parent, edge[1:]) in unique_chains:
				unique_chains.remove(chain)
				break

for node in nodelens:
	if node in long_nodes:
		print(node)
		continue
	chain = find(chain_parent, node)
	unique_chain = chain in unique_chains
	normalized_coverage = normalized_node_coverage[node]
	unique = False
	if unique_chain:
		chain_coverage = float(chain_coverage_sum[chain]) / float(chain_coverage_count[chain])
		chain_normalized_coverage = normalized_node_coverage[node] / chain_coverage
		if chain_normalized_coverage > 0.6 and chain_normalized_coverage < 1.4: unique = True
	nodelen = nodelens[node]
	if nodelen > 20000 and normalized_coverage > 0.9 and normalized_coverage < 1.1: unique = True
	if nodelen > 30000 and normalized_coverage > 0.8 and normalized_coverage < 1.2: unique = True
	if nodelen > 50000 and normalized_coverage > 0.7 and normalized_coverage < 1.3: unique = True
	if nodelen > 100000 and normalized_coverage > 0.6 and normalized_coverage < 1.4: unique = True
	if unique: print(node)
