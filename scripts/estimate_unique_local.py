#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
node_coverage_file = sys.argv[2]
alignment_file = sys.argv[3]
long_node_threshold = int(sys.argv[4])
solid_edge_threshold = int(sys.argv[5])
path_consistency_threshold = float(sys.argv[6])
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
	if len(edges[s]) == 1:
		if len(edges[revnode(getone(edges[s]))]) == 1:
			return getone(edges[s])
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
max_overlap = {}
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
			if fromnode not in max_overlap: max_overlap[fromnode] = int(parts[5][:-1])
			max_overlap[fromnode] = max(max_overlap[fromnode], int(parts[5][:-1]))
			if tonode not in max_overlap: max_overlap[tonode] = int(parts[5][:-1])
			max_overlap[tonode] = max(max_overlap[tonode], int(parts[5][:-1]))
			edges[fromnode].add(revnode(tonode))
			edges[tonode].add(revnode(fromnode))
			canon_edges.add(canontip(fromnode, tonode))

node_nonoverlap_lens = {}
for node in nodelens:
	overlaps = 0
	if ">" + node in max_overlap: overlaps += max_overlap[">" + node]
	if "<" + node in max_overlap: overlaps += max_overlap["<" + node]
	if overlaps >= nodelens[node]: overlaps = nodelens[node]-1
	node_nonoverlap_lens[node] = nodelens[node] - overlaps

for node in nodelens:
	parent[">" + node] = ">" + node
	rank[">" + node] = 1
	parent["<" + node] = "<" + node
	rank["<" + node] = 1
	if node not in long_nodes: merge(parent, rank, ">" + node, "<" + node)

for edge in canon_edges:
	merge(parent, rank, edge[0], edge[1])

cluster_edge_nodes = {}

for node in long_nodes:
	key = find(parent, ">" + node)
	if key not in cluster_edge_nodes: cluster_edge_nodes[key] = set()
	cluster_edge_nodes[key].add(">" + node)
	key = find(parent, "<" + node)
	if key not in cluster_edge_nodes: cluster_edge_nodes[key] = set()
	cluster_edge_nodes[key].add("<" + node)

chain_of_longnode = set()
for node in long_nodes:
	bubble_end = find_bubble_end(edges, ">" + node)
	while bubble_end is not None:
		next_node = bubble_end
		chain_of_longnode.add(next_node[1:])
		bubble_end = find_bubble_end(edges, next_node)
	bubble_end = find_bubble_end(edges, "<" + node)
	while bubble_end is not None:
		next_node = bubble_end
		chain_of_longnode.add(next_node[1:])
		bubble_end = find_bubble_end(edges, next_node)

node_coverage = {}
with open(node_coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node" and parts[1] == "coverage": continue
		assert parts[0] in nodelens
		assert parts[0] not in node_coverage
		node_coverage[parts[0]] = float(parts[1])

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
				if part_path[i] not in out_paths: out_paths[part_path[i]] = []
				if revnode(part_path[i]) not in out_paths: out_paths[revnode(part_path[i])] = []
				if i < len(part_path)-1: out_paths[part_path[i]].append(tuple(part_path[i+1:]))
				if i > 0: out_paths[revnode(part_path[i])].append(tuple(revnode(n) for n in part_path[0:i][::-1]))

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
			compare_coverage += float(node_coverage[end[1:]])
		compare_coverage /= float(len(compare_nodes))
	if compare_coverage < .1 * global_average_coverage or compare_coverage > 10 * global_average_coverage:
		sys.stderr.write("WARNING: nonsense local coverage of " + str(compare_coverage) + " for node " + node + ", reverting to global average " + str(global_average_coverage) + "\n")
		compare_coverage = global_average_coverage
	if compare_coverage > 0.01:
		normalized_node_coverage[node] = node_coverage[node] / compare_coverage
	else:
		sys.stderr.write('WARN: division by zero prevented for node ' + node + "\n")
		normalized_node_coverage[node] = 100.

roughly_average_coverage_nodes = set()
for n in nodelens:
	if n not in normalized_node_coverage: continue
	if normalized_node_coverage[n] < 0.6 and node_coverage[n] / global_average_coverage < 0.6: continue
	if normalized_node_coverage[n] > 1.4 and node_coverage[n] / global_average_coverage > 1.4: continue
	roughly_average_coverage_nodes.add(n)

path_consistent_nodes = set()
for n in nodelens:
	if n not in roughly_average_coverage_nodes: continue
	fw_paths = []
	bw_paths = []
	if ">" + n in out_paths: fw_paths = list(out_paths[">" + n])
	if "<" + n in out_paths: bw_paths = list(out_paths["<" + n])
	fw_paths.sort(key=lambda x: len(x))
	bw_paths.sort(key=lambda x: len(x))
	# if len(fw_paths) == 0: continue
	# if len(bw_paths) == 0: continue
	if len(fw_paths) + len(bw_paths) == 0: continue
	most_fw_consistent = 0
	for i in range(len(fw_paths)-1, -1, -1):
		consistents = 0
		inconsistents = 0
		for j in range(0, i):
			assert len(fw_paths[j]) <= len(fw_paths[i])
			if fw_paths[j] != fw_paths[i][0:len(fw_paths[j])]:
				inconsistents += 1
				continue
			consistents += 1
		if consistents > most_fw_consistent: most_fw_consistent = consistents
	if float(most_fw_consistent + len(bw_paths)) / float(len(fw_paths) + len(bw_paths)) < path_consistency_threshold: continue
	most_bw_consistent = 0
	for i in range(len(bw_paths)-1, -1, -1):
		consistents = 0
		inconsistents = 0
		for j in range(0, i):
			assert len(bw_paths[j]) <= len(bw_paths[i])
			if bw_paths[j] != bw_paths[i][0:len(bw_paths[j])]:
				inconsistents += 1
				continue
			consistents += 1
		if consistents > most_bw_consistent: most_bw_consistent = consistents
	if float(most_fw_consistent + most_bw_consistent) / float(len(fw_paths) + len(bw_paths)) < path_consistency_threshold: continue
	path_consistent_nodes.add(n)

path_unique_nodes = set()
for node in path_consistent_nodes:
	if node in roughly_average_coverage_nodes:
		path_unique_nodes.add(node)

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
	# if nodelens[node] < 20000: continue
	chain_coverage_count[chain] += nodelens[node]
	chain_coverage_sum[chain] += nodelens[node] * normalized_node_coverage[node]

copycount_2_chain_core_nodes = set()
for node in nodelens:
	chain = find(chain_parent, node)
	assert chain in chain_coverage_count
	assert chain in chain_coverage_sum
	if chain_coverage_count[chain] == 0: continue
	if chain_coverage_sum[chain] == 0: continue
	chain_coverage = float(chain_coverage_sum[chain]) / float(chain_coverage_count[chain])
	chain_length = float(chain_coverage_sum[chain])
	if chain_length > 75000 and chain_coverage > 1.5 and chain_coverage < 2.5: copycount_2_chain_core_nodes.add(node)

copycount_2_chain_unique_nodes = set()
copycount_2_chain_nonunique_nodes = set()
for node in copycount_2_chain_core_nodes:
	node_cov = normalized_node_coverage[node]
	if ">" + node in edges and len(edges[">" + node]) == 2:
		maybe_valid = True
		min_cov = 2
		max_cov = 0
		for edge in edges[">" + node]:
			if normalized_node_coverage[edge[1:]] < 0.5 or normalized_node_coverage[edge[1:]] > 1.5 or edge[1:] in copycount_2_chain_core_nodes: maybe_valid = False
			if len(edges[revnode(edge)]) != 1: maybe_valid = False
			min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
			max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
		if maybe_valid and max_cov < min_cov * 2 and max_cov >= min_cov:
			for edge in edges[">" + node]:
				copycount_2_chain_unique_nodes.add(edge[1:])
		if maybe_valid and max_cov > min_cov * 2:
			for edge in edges[">" + node]:
				copycount_2_chain_nonunique_nodes.add(edge[1:])
	if "<" + node in edges and len(edges["<" + node]) == 2:
		maybe_valid = True
		min_cov = 2
		max_cov = 0
		for edge in edges["<" + node]:
			if normalized_node_coverage[edge[1:]] < 0.5 or normalized_node_coverage[edge[1:]] > 1.5 or edge[1:] in copycount_2_chain_core_nodes: maybe_valid = False
			if len(edges[revnode(edge)]) != 1: maybe_valid = False
			min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
			max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
		if maybe_valid and max_cov < min_cov * 2 and max_cov >= min_cov:
			for edge in edges["<" + node]:
				copycount_2_chain_unique_nodes.add(edge[1:])
		if maybe_valid and max_cov > min_cov * 2:
			for edge in edges["<" + node]:
				copycount_2_chain_nonunique_nodes.add(edge[1:])
	if ">" + node in edges and len(edges[">" + node]) == 3:
		maybe_valid = True
		further_ahead_nodes = {}
		maybe_addable = set()
		min_cov = 2
		max_cov = 0
		valid_direct = 0
		valid_indirect = 0
		maybe_disablable = set()
		for edge in edges[">" + node]:
			if edge in edges and len(edges[edge]) == 1:
				for edge2 in edges[edge]:
					if len(edges[revnode(edge2)]) == 2 and edge2[1:] not in copycount_2_chain_core_nodes:
						if edge2 not in further_ahead_nodes: further_ahead_nodes[edge2] = 0
						further_ahead_nodes[edge2] += 1
		if len(further_ahead_nodes) <= 2:
			max_further_ahead = 0
			min_further_ahead = 3
			for edge2 in further_ahead_nodes:
				max_further_ahead = max(max_further_ahead, further_ahead_nodes[edge2])
				min_further_ahead = min(min_further_ahead, further_ahead_nodes[edge2])
			if (min_further_ahead == 1 and max_further_ahead == 2) or (max_further_ahead == 2 and len(further_ahead_nodes) == 1):
				for edge in edges[">" + node]:
					skip_this = False
					if edge in edges and len(edges[edge]) == 1:
						for edge2 in edges[edge]:
							if edge2 in further_ahead_nodes and further_ahead_nodes[edge2] == 2:
								skip_this = True
					if skip_this:
						maybe_disablable.add(edge[1:])
					else:
						if normalized_node_coverage[edge[1:]] > 0.5 and normalized_node_coverage[edge[1:]] < 1.5 and edge[1:] not in copycount_2_chain_core_nodes:
							maybe_addable.add(edge[1:])
							min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
							max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
				for edge in further_ahead_nodes:
					if further_ahead_nodes[edge] == 2:
						if normalized_node_coverage[edge[1:]] > 0.5 and normalized_node_coverage[edge[1:]] < 1.5 and edge[1:] not in copycount_2_chain_core_nodes:
							maybe_addable.add(edge[1:])
							min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
							max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
			else:
				maybe_valid = False
		else:
			maybe_valid = False
		if maybe_valid and max_cov < min_cov * 2 and max_cov >= min_cov:
			for n in maybe_addable: copycount_2_chain_unique_nodes.add(n)
			for n in maybe_disablable: copycount_2_chain_nonunique_nodes.add(n)
	if "<" + node in edges and len(edges["<" + node]) == 3:
		maybe_valid = True
		further_ahead_nodes = {}
		maybe_addable = set()
		min_cov = 2
		max_cov = 0
		valid_direct = 0
		valid_indirect = 0
		maybe_disablable = set()
		for edge in edges["<" + node]:
			if edge in edges and len(edges[edge]) == 1:
				for edge2 in edges[edge]:
					if len(edges[revnode(edge2)]) == 2 and edge2[1:] not in copycount_2_chain_core_nodes:
						if edge2 not in further_ahead_nodes: further_ahead_nodes[edge2] = 0
						further_ahead_nodes[edge2] += 1
		if len(further_ahead_nodes) <= 2:
			max_further_ahead = 0
			min_further_ahead = 3
			for edge2 in further_ahead_nodes:
				max_further_ahead = max(max_further_ahead, further_ahead_nodes[edge2])
				min_further_ahead = min(min_further_ahead, further_ahead_nodes[edge2])
			if (min_further_ahead == 1 and max_further_ahead == 2) or (max_further_ahead == 2 and len(further_ahead_nodes) == 1):
				for edge in edges["<" + node]:
					skip_this = False
					if edge in edges and len(edges[edge]) == 1:
						for edge2 in edges[edge]:
							if edge2 in further_ahead_nodes and further_ahead_nodes[edge2] == 2:
								skip_this = True
					if skip_this:
						maybe_disablable.add(edge[1:])
					else:
						if normalized_node_coverage[edge[1:]] > 0.5 and normalized_node_coverage[edge[1:]] < 1.5 and edge[1:] not in copycount_2_chain_core_nodes:
							maybe_addable.add(edge[1:])
							min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
							max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
				for edge in further_ahead_nodes:
					if further_ahead_nodes[edge] == 2:
						if normalized_node_coverage[edge[1:]] > 0.5 and normalized_node_coverage[edge[1:]] < 1.5 and edge[1:] not in copycount_2_chain_core_nodes:
							maybe_addable.add(edge[1:])
							min_cov = min(min_cov, normalized_node_coverage[edge[1:]])
							max_cov = max(max_cov, normalized_node_coverage[edge[1:]])
			else:
				maybe_valid = False
		else:
			maybe_valid = False
		if maybe_valid and max_cov < min_cov * 2 and max_cov >= min_cov:
			for n in maybe_addable: copycount_2_chain_unique_nodes.add(n)
			for n in maybe_disablable: copycount_2_chain_nonunique_nodes.add(n)


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

chain_unique_nodes = set()
length_unique_nodes = set()
for node in nodelens:
	if node in long_nodes:
		continue
	chain = find(chain_parent, node)
	unique_chain = chain in unique_chains
	normalized_coverage = normalized_node_coverage[node]
	unique = False
	if unique_chain:
		chain_coverage = float(chain_coverage_sum[chain]) / float(chain_coverage_count[chain])
		chain_normalized_coverage = normalized_node_coverage[node] / chain_coverage
		if chain_normalized_coverage > 0.6 and chain_normalized_coverage < 1.4: chain_unique_nodes.add(node)
	nodelen = node_nonoverlap_lens[node]
	if nodelen > 20000 and normalized_coverage > 0.9 and normalized_coverage < 1.1: length_unique_nodes.add(node)
	if nodelen > 30000 and normalized_coverage > 0.8 and normalized_coverage < 1.2: length_unique_nodes.add(node)
	if nodelen > 50000 and normalized_coverage > 0.7 and normalized_coverage < 1.3: length_unique_nodes.add(node)
	if nodelen > 100000 and normalized_coverage > 0.6 and normalized_coverage < 1.4: length_unique_nodes.add(node)

uniques = set()
for node in long_nodes: uniques.add(node)
for node in path_unique_nodes: uniques.add(node)
for node in chain_of_longnode: uniques.add(node)
for node in copycount_2_chain_unique_nodes: uniques.add(node)
for node in chain_unique_nodes: uniques.add(node)
for node in length_unique_nodes: uniques.add(node)
for node in copycount_2_chain_nonunique_nodes: 
	if node in uniques: uniques.remove(node)
for node in copycount_2_chain_core_nodes: 
	if node in uniques: uniques.remove(node)

for node in uniques:
	print(node)
