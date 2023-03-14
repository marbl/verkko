#!/usr/bin/env python

import sys
import heapq

input_gfa = sys.argv[1]
out_path_file = sys.argv[2]
node_coverage_file = sys.argv[3]
resolve_namemapping_file = sys.argv[4]
max_resolve_length = int(sys.argv[5])
min_allowed_coverage = float(sys.argv[6])
resolve_steps = [int(n) for n in sys.argv[7:]]
# gaf from stdin
# gfa to stdout

next_unitig_num = 1

def iterate_deterministic(l):
	tmp = list(l)
	tmp.sort()
	for x in tmp:
		yield x

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
	return "".join(comp[c] for c in s[::-1])

def getone(s):
	assert len(s) == 1
	for n in s:
		return n

def revnode(n):
	assert len(n) >= 2
	assert n[0] == ">" or n[0] == "<"
	return (">" if n[0] == "<" else "<") + n[1:]

def canon(node1, node2):
	fwstr = node1 + node2
	bwstr = revnode(node2) + revnode(node1)
	if bwstr < fwstr: return (revnode(node2), revnode(node1))
	return (node1, node2)

def read_graph(filename):
	base_seqs = {}
	node_seqs = {}
	edges = {}
	edge_overlaps = {}
	with open(filename) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == 'S':
				base_seqs[parts[1]] = parts[2]
				node_seqs[parts[1]] = ([">" + parts[1]], 0, 0)
			if parts[0] == 'L':
				fromnode = (">" if parts[2] == "+" else "<") + parts[1]
				tonode = (">" if parts[4] == "+" else "<") + parts[3]
				if fromnode not in edges: edges[fromnode] = set()
				edges[fromnode].add(tonode)
				if revnode(tonode) not in edges: edges[revnode(tonode)] = set()
				edges[revnode(tonode)].add(revnode(fromnode))
				key = canon(fromnode, tonode)
				edge_overlaps[key] = int(parts[5][:-1])
	return (base_seqs, node_seqs, edges, edge_overlaps)

def read_graph_only_bases(filename):
	base_seqs = {}
	with open(filename) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == 'S':
				base_seqs[parts[1]] = parts[2]
	return base_seqs

def replace_path_nodes(resolvable, paths_crossing, new_edgenodes, maybe_resolvable):
	found_identities = set()
	remove_paths = []
	add_paths = []
	for node in iterate_deterministic(resolvable):
		for path in iterate_paths(paths_crossing, node):
			if id(path) in found_identities: continue
			found_identities.add(id(path))
			remove_paths.append(path)
			if len(path) < 2:
				continue
			new_path = []
			for j in range(0, len(path)):
				if path[j][1:] not in resolvable:
					new_path.append(path[j])
				if j > 0 and (revnode(path[j]), revnode(path[j-1])) in new_edgenodes:
					new_path.append("<" + new_edgenodes[(revnode(path[j]), revnode(path[j-1]))])
				if j < len(path)-1 and (path[j], path[j+1]) in new_edgenodes:
					new_path.append(">" + new_edgenodes[(path[j], path[j+1])])
			add_paths.append(new_path)
	for path in remove_paths:
		remove_path(paths_crossing, path)
	for path in add_paths:
		for node in path: assert node[1:] not in resolvable
		if len(path) > 1: add_path(paths_crossing, path)
	for node in resolvable: assert len(paths_crossing[node]) == 0

def replace_path_node(paths_crossing, node, left_added, right_added):
	remove_paths = []
	add_paths = []
	for path in iterate_paths(paths_crossing, node):
		remove_paths.append(path)
		if len(path) < 2:
			continue
		new_path = []
		for j in range(0, len(path)):
			if path[j][1:] != node:
				new_path.append(path[j])
				continue
			assert path[j] == ">" + node or path[j] == "<" + node
			if path[j] == ">" + node:
				if j > 0:
					assert path[j-1] in left_added
					new_path.append(">" + left_added[path[j-1]])
				if j < len(path)-1:
					assert path[j+1] in right_added
					new_path.append(">" + right_added[path[j+1]])
			else:
				assert path[j] == "<" + node
				if j > 0:
					assert revnode(path[j-1]) in right_added
					new_path.append("<" + right_added[revnode(path[j-1])])
				if j < len(path)-1:
					assert revnode(path[j+1]) in left_added
					new_path.append("<" + left_added[revnode(path[j+1])])
		add_paths.append(new_path)
	for path in remove_paths:
		remove_path(paths_crossing, path)
	assert len(paths_crossing[node]) == 0
	for path in add_paths:
		assert node not in path
		if len(path) > 1: add_path(paths_crossing, path)
	assert len(paths_crossing[node]) == 0

def remove_graph_node(node, node_seqs, edges):
	assert node in node_seqs
	del node_seqs[node]
	remove_edges = []
	if ">" + node in edges:
		for edge in iterate_deterministic(edges[">" + node]):
			if edge[1:] == node: continue
			assert revnode(edge) in edges
			assert "<" + node in edges[revnode(edge)]
			remove_edges.append((revnode(edge), "<" + node))
		del edges[">" + node]
	if "<" + node in edges:
		for edge in iterate_deterministic(edges["<" + node]):
			if edge[1:] == node: continue
			assert revnode(edge) in edges
			assert ">" + node in edges[revnode(edge)]
			remove_edges.append((revnode(edge), ">" + node))
		del edges["<" + node]
	for edge in remove_edges:
		if edge[0] not in edges: continue
		edges[edge[0]].remove(edge[1])

def print_relevant_paths(node, paths_crossing):
	for path in iterate_paths(paths_crossing, node):
		sys.stderr.write(str(path) + "\n")

def get_valid_triplets(node, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes, node_seqs):
	if ">" + node not in edges and "<" not in edges: return []
	if ">" + node not in edges: edges[">" + node] = set()
	if "<" + node not in edges: edges["<" + node] = set()
	if len(edges[">" + node]) <= 1 and len(edges["<" + node]) <= 1: return []
	if len(edges[">" + node]) == 1 and getone(edges[">" + node]) == "<" + node: return []
	if len(edges["<" + node]) == 1 and getone(edges["<" + node]) == ">" + node: return []
	triplets = {}
	covered_in_neighbors = {}
	covered_out_neighbors = {}
	if node not in paths_crossing: return []
	for p in iterate_paths(paths_crossing, node):
		if len(p) == 1: continue
		for i in range(1, len(p)):
			if p[i] == ">" + node:
				if p[i-1] not in covered_in_neighbors: covered_in_neighbors[p[i-1]] = 0
				covered_in_neighbors[p[i-1]] += 1
			if p[i] == "<" + node:
				if revnode(p[i-1]) not in covered_out_neighbors: covered_out_neighbors[revnode(p[i-1])] = 0
				covered_out_neighbors[revnode(p[i-1])] += 1
			if p[i-1] == ">" + node:
				if p[i] not in covered_out_neighbors: covered_out_neighbors[p[i]] = 0
				covered_out_neighbors[p[i]] += 1
			if p[i-1] == "<" + node:
				if revnode(p[i]) not in covered_in_neighbors: covered_in_neighbors[revnode(p[i])] = 0
				covered_in_neighbors[revnode(p[i])] += 1
		for i in range(0, len(p)):
			if p[i][1:] != node: continue
			if i == 0 and len(edges[revnode(p[0])]) != 0: continue
			if i == len(p)-1 and len(edges[p[i]]) != 0: continue
			triplet = p[i-1:i+2]
			if i == 0: triplet = [None, p[0], p[1]]
			if i == len(p)-1: triplet = [p[i-1], p[i], None]
			assert triplet[1] == ">" + node or triplet[1] == "<" + node
			if triplet[1] == "<" + node:
				if i != 0 and i != len(p)-1:
					triplet = [revnode(n) for n in triplet[::-1]]
				elif i == 0:
					triplet = (revnode(p[1]), revnode(p[0]), None)
				else:
					assert i == len(p)-1
					triplet = (None, revnode(p[i]), revnode(p[i-1]))
			assert triplet[1] == ">" + node
			triplet = tuple(triplet)
			if triplet not in triplets: triplets[triplet] = 0
			triplets[triplet] += 1
	if len(triplets) == 0: return []
	triplet_covered_in_neighbors = set()
	triplet_covered_out_neighbors = set()
	solid_triplets = set()
	for triplet in triplets:
		(fromnode, middle, tonode) = triplet
		cov = triplets[triplet]
		assert middle == ">" + node
		if cov < min_edge_support: continue
		solid_triplets.add(triplet)
	for triplet in solid_triplets:
		(fromnode, middle, tonode) = triplet
		if fromnode is not None: triplet_covered_in_neighbors.add(fromnode)
		if tonode is not None: triplet_covered_out_neighbors.add(tonode)
	for edge in edges[">" + node]:
		cov = 0
		if edge in covered_out_neighbors: cov = covered_out_neighbors[edge]
		removable = True
		for partnode in node_seqs[edge[1:]][0]:
			if partnode[1:] not in removable_nodes: removable = False
		if cov < min_coverage and removable: continue
		if edge not in triplet_covered_out_neighbors: return []
	for edge in edges["<" + node]:
		cov = 0
		if revnode(edge) in covered_in_neighbors: cov = covered_in_neighbors[revnode(edge)]
		removable = True
		for partnode in node_seqs[edge[1:]][0]:
			if partnode[1:] not in removable_nodes: removable = False
		if cov < min_coverage and removable: continue
		if revnode(edge) not in triplet_covered_in_neighbors: return []
	allowed_triplets = list(iterate_deterministic(solid_triplets))
	return allowed_triplets

def resolve_hairpins(nodelength, nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, min_edge_support, min_coverage, removable_nodes):
	hairpins = set()
	for node in iterate_deterministic(nodes):
		if ">" + node not in edges: continue
		if "<" + node not in edges: continue
		if len(edges[">" + node]) == 1 and getone(edges[">" + node]) == "<" + node:
			hairpins.add(">" + node)
		if len(edges["<" + node]) == 1 and getone(edges["<" + node]) == ">" + node:
			hairpins.add("<" + node)
	remove_paths = []
	add_paths = []
	new_node_names = set()
	resolved = set()
	for node in hairpins:
		assert revnode(node) not in hairpins # double hairpin resolution hard to implement, so just hope it never happens
		resolutions = {}
		for path in iterate_paths(paths_crossing, node[1:]):
			if len(path) < 4: continue
			for i in range(1, len(path)-2):
				if path[i] == node and path[i+1] == revnode(node):
					key = canon(path[i-1], path[i+2])
					if key not in resolutions: resolutions[key] = 0
					resolutions[key] += 1
				if path[i] == revnode(node) and path[i+1] == node:
					assert False # this should never happen?
		covered_edges = set()
		solid_resolutions = set()
		for resolution in resolutions:
			if resolutions[resolution] < min_edge_support: continue
			solid_resolutions.add(resolution)
			covered_edges.add(revnode(resolution[0]))
			covered_edges.add(resolution[1])
		if len(covered_edges) < len(edges[revnode(node)]): continue
		sys.stderr.write("resolve hairpin " + node + "\n")
		resolved.add(node[1:])
		nextnum = 0
		resolution_number = {}
		for key in solid_resolutions:
			fwname = node[1:] + "hairpin" + str(nextnum) + "fw"
			bwname = node[1:] + "hairpin" + str(nextnum) + "bw"
			new_node_names.add(fwname)
			new_node_names.add(bwname)
			resolution_number[key] = nextnum
			nextnum += 1
			if node[0] == ">":
				node_seqs[fwname] = node_seqs[node[1:]]
				node_seqs[bwname] = ([revnode(n) for n in node_seqs[node[1:]][0][::-1]], node_seqs[node[1:]][2], node_seqs[node[1:]][1])
			else:
				node_seqs[fwname] = ([revnode(n) for n in node_seqs[node[1:]][0][::-1]], node_seqs[node[1:]][2], node_seqs[node[1:]][1])
				node_seqs[bwname] = node_seqs[node[1:]]
			edges["<" + fwname] = set()
			edges[">" + fwname] = set()
			edges["<" + bwname] = set()
			edges[">" + bwname] = set()
			edges[">" + fwname].add(">" + bwname)
			edges["<" + bwname].add("<" + fwname)
			edges["<" + fwname].add(revnode(key[0]))
			edges[">" + bwname].add(key[1])
			edges[key[0]].add(">" + fwname)
			edges[revnode(key[1])].add("<" + bwname)
			edge_overlaps[canon(">" + fwname, ">" + bwname)] = edge_overlaps[canon(node, revnode(node))]
			edge_overlaps[canon(key[0], ">" + fwname)] = edge_overlaps[canon(key[0], node)]
			edge_overlaps[canon(">" + bwname, key[1])] = edge_overlaps[canon(revnode(node), key[1])]
		remove_graph_node(node[1:], node_seqs, edges)
		for path in iterate_paths(paths_crossing, node[1:]):
			remove_paths.append(path)
			if len(path) < 4: continue
			add_this = []
			add_this.append(path[0])
			for i in range(1, len(path)-2):
				if path[i] == node and path[i+1] == revnode(node):
					key = (path[i-1], path[i+2])
					canonkey = canon(key[0], key[1])
					if canonkey not in resolution_number:
						add_paths.append(add_this)
						add_this = []
						continue
					if key == canonkey:
						add_this.append(">" + node[1:] + "hairpin" + str(resolution_number[canonkey]) + "fw")
						add_this.append(">" + node[1:] + "hairpin" + str(resolution_number[canonkey]) + "bw")
					else:
						add_this.append("<" + node[1:] + "hairpin" + str(resolution_number[canonkey]) + "bw")
						add_this.append("<" + node[1:] + "hairpin" + str(resolution_number[canonkey]) + "fw")
				elif path[i] == revnode(node) and path[i+1] == node:
					assert False # this should never happen?
				elif path[i][1:] == node[1:]:
					continue
				else:
					add_this.append(path[i])
			if path[-2][1:] != node[1:]: add_this.append(path[-2])
			if path[-1][1:] != node[1:]: add_this.append(path[-1])
			add_paths.append(add_this)
	for path in remove_paths:
		remove_path(paths_crossing, path)
	for path in add_paths:
		add_path(paths_crossing, path)
	return (new_node_names, resolved)

def resolve_nodes(nodelength, nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, min_edge_support, min_coverage, removable_nodes):
	resolvable = set()
	triplets = []
	for node in iterate_deterministic(nodes):
		triplets_here = get_valid_triplets(node, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes, node_seqs)
		if len(triplets_here) == 0:
			maybe_resolvable.remove(node)
			# print("no triplets at " + node)
		else:
			resolvable.add(node)
			triplets += triplets_here
	fixed_any = True
	while fixed_any:
		fixed_any = False
		for triplet in triplets:
			if triplet[1][1:] not in maybe_resolvable: continue
			if triplet[0] is not None and not (triplet[0][1:] not in maybe_resolvable and get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[0][1:]) - edge_overlaps[canon(triplet[0], triplet[1])] == 1): continue
			if triplet[2] is not None and not (triplet[2][1:] not in maybe_resolvable and get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[2][1:]) - edge_overlaps[canon(triplet[1], triplet[2])] == 1): continue
			maybe_resolvable.remove(triplet[1][1:])
			# print("borders unresolvable at " + triplet[1][1:])
			fixed_any = True
	resolvable = set()
	triplets = []
	longest_extension_per_node = {}
	for node in iterate_deterministic(nodes):
		if node not in maybe_resolvable: continue
		triplets_here = get_valid_triplets(node, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes, node_seqs)
		if len(triplets_here) == 0:
			maybe_resolvable.remove(node)
		else:
			resolvable.add(node)
			triplets += triplets_here
			# use max step size 1000 to avoid strange multiple megabase bubbles
			if ">" + node not in longest_extension_per_node: longest_extension_per_node[">" + node] = 1000
			if "<" + node not in longest_extension_per_node: longest_extension_per_node["<" + node] = 1000
			for triplet in triplets_here:
				if triplet[0] is not None:
					left_node_length = get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[0][1:]) - edge_overlaps[canon(triplet[0], triplet[1])]
					assert left_node_length >= 1
					if not (triplet[0][1:] not in maybe_resolvable and left_node_length == 1):
						if left_node_length > 1: left_node_length -= 1
						if left_node_length < longest_extension_per_node["<" + node]: longest_extension_per_node["<" + node] = left_node_length
				if triplet[2] is not None:
					right_node_length = get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[2][1:]) - edge_overlaps[canon(triplet[1], triplet[2])]
					assert right_node_length >= 1
					if not (triplet[2][1:] not in maybe_resolvable and right_node_length == 1):
						if right_node_length > 1: right_node_length -= 1
						if right_node_length < longest_extension_per_node[">" + node]: longest_extension_per_node[">" + node] = right_node_length
			assert longest_extension_per_node[">" + node] > 0
			assert longest_extension_per_node["<" + node] > 0
	new_edgenodes = {}
	for triplet in triplets:
		assert triplet[1][1:] in maybe_resolvable
		assert triplet[1][0] == ">"
		extend_amount = longest_extension_per_node["<" + triplet[1][1:]]
		if triplet[0] is not None:
			key = (revnode(triplet[1]), revnode(triplet[0]))
			if key[1][1:] in maybe_resolvable or get_unitig_len(node_lens, edge_overlaps, node_seqs, key[1][1:]) - edge_overlaps[canon(key[0], key[1])] > 1:
				nodename = "edge_" + str(key[0][1:]) + ("fw" if key[0][0] == ">" else "bw") + "_" + str(key[1][1:]) + ("fw" if key[1][0] == ">" else "bw")
				new_edgenodes[key] = nodename
				assert key[0][0] == "<"
				if key[0][0] == ">":
					node_seqs[nodename] = (node_seqs[key[0][1:]][0], node_seqs[key[0][1:]][1], node_seqs[key[0][1:]][2] - extend_amount)
				else:
					node_seqs[nodename] = ([revnode(n) for n in node_seqs[key[0][1:]][0][::-1]], node_seqs[key[0][1:]][2], node_seqs[key[0][1:]][1] - extend_amount)
				if node_seqs[nodename][2] < 0:
					last_node_seq = node_seqs[key[0][1:]][0]
					new_node_seq = node_seqs[key[1][1:]][0]
					last_end_clip = node_seqs[key[0][1:]][2]
					new_start_clip = node_seqs[key[1][1:]][1]
					last_start_clip = node_seqs[key[0][1:]][1]
					new_end_clip = node_seqs[key[1][1:]][2]
					if key[0][0] == "<":
						last_node_seq = [revnode(n) for n in last_node_seq[::-1]]
						last_end_clip = node_seqs[key[0][1:]][1]
						last_start_clip = node_seqs[key[0][1:]][2]
					if key[1][0] == "<":
						new_node_seq = [revnode(n) for n in new_node_seq[::-1]]
						new_start_clip = node_seqs[key[1][1:]][2]
						new_end_clip = node_seqs[key[1][1:]][1]
					overlap = edge_overlaps[canon(key[0], key[1])]
					found = False
					for i in range(0, len(last_node_seq)):
						if node_lens[last_node_seq[i][1:]] <= new_start_clip: continue
						if get_unitig_len_path(node_lens, edge_overlaps, last_node_seq[i:], new_start_clip, last_end_clip) == overlap:
							assert last_node_seq[i:] == new_node_seq[:len(last_node_seq)-i]
							assert i > 0 or last_start_clip < new_start_clip
							result_node_seq = last_node_seq + new_node_seq[len(last_node_seq)-i:]
							found = True
							break
					if not found:
						assert edge_overlaps[canon(last_node_seq[-1], new_node_seq[0])] == overlap
						result_node_seq = last_node_seq + new_node_seq
					assert result_node_seq[:len(last_node_seq)] == last_node_seq
					new_length = get_unitig_len_path(node_lens, edge_overlaps, result_node_seq, last_start_clip, new_end_clip)
					old_length = get_unitig_len(node_lens, edge_overlaps, node_seqs, key[0][1:])
					assert new_length >= old_length + extend_amount
					result_end_clip = new_end_clip + (new_length - (old_length + extend_amount))
					while result_end_clip >= node_lens[result_node_seq[-1][1:]]:
						assert len(result_node_seq) >= 2
						result_end_clip -= node_lens[result_node_seq[-1][1:]] - edge_overlaps[canon(result_node_seq[-2], result_node_seq[-1])]
						result_node_seq = result_node_seq[:-1]
					assert result_end_clip >= 0
					node_seqs[nodename] = (result_node_seq, last_start_clip, result_end_clip)
				if key[1][1:] in resolvable:
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, nodename) == nodelength+extend_amount
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, key[0][1:]) == nodelength
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, key[1][1:]) == nodelength
		extend_amount = longest_extension_per_node[">" + triplet[1][1:]]
		if triplet[2] is not None:
			key = (triplet[1], triplet[2])
			if key[1][1:] in maybe_resolvable or get_unitig_len(node_lens, edge_overlaps, node_seqs, key[1][1:]) - edge_overlaps[canon(key[0], key[1])] > 1:
				nodename = "edge_" + str(key[0][1:]) + ("fw" if key[0][0] == ">" else "bw") + "_" + str(key[1][1:]) + ("fw" if key[1][0] == ">" else "bw")
				new_edgenodes[key] = nodename
				assert key[0][0] == ">"
				if key[0][0] == ">":
					node_seqs[nodename] = (node_seqs[key[0][1:]][0], node_seqs[key[0][1:]][1], node_seqs[key[0][1:]][2] - extend_amount)
				else:
					node_seqs[nodename] = ([revnode(n) for n in node_seqs[key[0][1:]][0][::-1]], node_seqs[key[0][1:]][2], node_seqs[key[0][1:]][1] - extend_amount)
				if node_seqs[nodename][2] < 0:
					last_node_seq = node_seqs[key[0][1:]][0]
					new_node_seq = node_seqs[key[1][1:]][0]
					last_end_clip = node_seqs[key[0][1:]][2]
					new_start_clip = node_seqs[key[1][1:]][1]
					last_start_clip = node_seqs[key[0][1:]][1]
					new_end_clip = node_seqs[key[1][1:]][2]
					if key[0][0] == "<":
						last_node_seq = [revnode(n) for n in last_node_seq[::-1]]
						last_end_clip = node_seqs[key[0][1:]][1]
						last_start_clip = node_seqs[key[0][1:]][2]
					if key[1][0] == "<":
						new_node_seq = [revnode(n) for n in new_node_seq[::-1]]
						new_start_clip = node_seqs[key[1][1:]][2]
						new_end_clip = node_seqs[key[1][1:]][1]
					overlap = edge_overlaps[canon(key[0], key[1])]
					found = False
					for i in range(0, len(last_node_seq)):
						if node_lens[last_node_seq[i][1:]] <= new_start_clip: continue
						if get_unitig_len_path(node_lens, edge_overlaps, last_node_seq[i:], new_start_clip, last_end_clip) == overlap:
							assert last_node_seq[i:] == new_node_seq[:len(last_node_seq)-i]
							assert i > 0 or last_start_clip < new_start_clip
							result_node_seq = last_node_seq + new_node_seq[len(last_node_seq)-i:]
							found = True
							break
					if not found:
						assert edge_overlaps[canon(last_node_seq[-1], new_node_seq[0])] == overlap
						result_node_seq = last_node_seq + new_node_seq
					assert result_node_seq[:len(last_node_seq)] == last_node_seq
					new_length = get_unitig_len_path(node_lens, edge_overlaps, result_node_seq, last_start_clip, new_end_clip)
					old_length = get_unitig_len(node_lens, edge_overlaps, node_seqs, key[0][1:])
					assert new_length >= old_length + extend_amount
					result_end_clip = new_end_clip + (new_length - (old_length + extend_amount))
					while result_end_clip >= node_lens[result_node_seq[-1][1:]]:
						assert len(result_node_seq) >= 2
						result_end_clip -= node_lens[result_node_seq[-1][1:]] - edge_overlaps[canon(result_node_seq[-2], result_node_seq[-1])]
						result_node_seq = result_node_seq[:-1]
					assert result_end_clip >= 0
					node_seqs[nodename] = (result_node_seq, last_start_clip, result_end_clip)
				if key[1][1:] in resolvable:
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, nodename) == nodelength+extend_amount
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, key[0][1:]) == nodelength
					assert get_unitig_len(node_lens, edge_overlaps, node_seqs, key[1][1:]) == nodelength
	new_node_names = set()
	for key in new_edgenodes:
		nodename = new_edgenodes[key]
		# if not get_unitig_len(node_lens, edge_overlaps, node_seqs, nodename) == nodelength + longest_extension_per_node[key[0][1:]]:
		# 	print(key)
		# 	print(node_seqs[key[0][1:]])
		# 	print(get_unitig_len(node_lens, edge_overlaps, node_seqs, key[0][1:]))
		# 	print(node_seqs[key[1][1:]])
		# 	print(get_unitig_len(node_lens, edge_overlaps, node_seqs, key[1][1:]))
		# 	print(node_seqs[nodename])
		# 	print(nodelength)
		# 	print(longest_extension_per_node[key[0][1:]])
		# 	print(get_unitig_len(node_lens, edge_overlaps, node_seqs, nodename))
		# assert get_unitig_len(node_lens, edge_overlaps, node_seqs, nodename) == nodelength + longest_extension_per_node[key[0][1:]]
		assert node_seqs[nodename][1] >= 0
		assert node_seqs[nodename][2] >= 0
		edges[">" + nodename] = set()
		edges["<" + nodename] = set()
		for i in range(1, len(node_seqs[nodename][0])):
			assert canon(node_seqs[nodename][0][i-1], node_seqs[nodename][0][i]) in edge_overlaps
		new_node_names.add(new_edgenodes[key])
	for triplet in triplets:
		overlap = nodelength
		assert triplet[1][1:] in resolvable
		has_left_key = False
		has_right_key = False
		extend_amount = longest_extension_per_node["<" + triplet[1][1:]]
		if triplet[0] is not None:
			if triplet[0][1:] not in maybe_resolvable and get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[0][1:]) - edge_overlaps[canon(triplet[0], triplet[1])] == 1:
				assert (revnode(triplet[1]), revnode(triplet[0])) not in new_edgenodes
				assert get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[0][1:]) <= nodelength
				left_key = triplet[0]
				overlap = edge_overlaps[canon(triplet[0], triplet[1])]
				has_left_key = False
			else:
				has_left_key = True
				assert (revnode(triplet[1]), revnode(triplet[0])) in new_edgenodes
				left_key = "<" + new_edgenodes[(revnode(triplet[1]), revnode(triplet[0]))]
				add_overlap = extend_amount
				if triplet[0][1:] not in resolvable:
					lefter_key = triplet[0]
				else:
					assert (triplet[0], triplet[1]) in new_edgenodes
					lefter_key = ">" + new_edgenodes[(triplet[0], triplet[1])]
					add_overlap = extend_amount + longest_extension_per_node[triplet[0]]
				edges[lefter_key].add(left_key)
				edges[revnode(left_key)].add(revnode(lefter_key))
				key = canon(lefter_key, left_key)
				edge_overlaps[key] = edge_overlaps[canon(triplet[0], triplet[1])]+add_overlap
		extend_amount = longest_extension_per_node[">" + triplet[1][1:]]
		if triplet[2] is not None:
			if triplet[2][1:] not in maybe_resolvable and get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[2][1:]) - edge_overlaps[canon(triplet[1], triplet[2])] == 1:
				assert (triplet[1], triplet[2]) not in new_edgenodes
				assert get_unitig_len(node_lens, edge_overlaps, node_seqs, triplet[2][1:]) <= nodelength
				right_key = triplet[2]
				overlap = edge_overlaps[canon(triplet[1], triplet[2])]
				has_right_key = False
			else:
				has_right_key = True
				assert (triplet[1], triplet[2]) in new_edgenodes
				right_key = ">" + new_edgenodes[(triplet[1], triplet[2])]
				add_overlap = extend_amount
				if triplet[2][1:] not in resolvable:
					righter_key = triplet[2]
				else:
					assert (revnode(triplet[2]), revnode(triplet[1])) in new_edgenodes
					righter_key = "<" + new_edgenodes[(revnode(triplet[2]), revnode(triplet[1]))]
					add_overlap = extend_amount + longest_extension_per_node[revnode(triplet[2])]
				edges[right_key].add(righter_key)
				edges[revnode(righter_key)].add(revnode(right_key))
				key = canon(right_key, righter_key)
				edge_overlaps[key] = edge_overlaps[canon(triplet[1], triplet[2])]+add_overlap
		if triplet[0] is None or triplet[2] is None: continue
		assert has_left_key or has_right_key
		edges[left_key].add(right_key)
		edges[revnode(right_key)].add(revnode(left_key))
		key = canon(left_key, right_key)
		assert key not in edge_overlaps
		edge_overlaps[key] = overlap
	replace_path_nodes(resolvable, paths_crossing, new_edgenodes, maybe_resolvable)
	for node in iterate_deterministic(resolvable): remove_graph_node(node, node_seqs, edges)
	self_edge_coverage = {}
	new_paths = []
	found_identities = set()
	for n in iterate_deterministic(new_edgenodes):
		for path in iterate_paths(paths_crossing, new_edgenodes[n]):
			if id(path) in found_identities: continue
			found_identities.add(id(path))
			new_paths.append(path)
	split_add_paths = []
	remove_add_paths = []
	for path in new_paths:
		last_split = 0
		for i in range(1, len(path)):
			if path[i-1] not in edges or path[i] not in edges[path[i-1]]:
				split_add_paths.append(path[last_split:i])
				last_split = i
			else:
				assert canon(path[i-1], path[i]) in edge_overlaps
		if last_split != 0:
			split_add_paths.append(path[last_split:])
			remove_add_paths.append(path)
	for path in remove_add_paths:
		remove_path(paths_crossing, path)
	for path in split_add_paths:
		add_path(paths_crossing, path)
	return (new_node_names, resolvable)

def find(parent, key):
	while parent[key] != parent[parent[key]]: parent[key] = parent[parent[key]]
	return parent[key]

def merge(parent, left, right):
	left = find(parent, left)
	right = find(parent, right)
	parent[right] = left

def add_safe_unitig(node, edges, safe_nodes, safe_edges):
	if node[1:] in safe_nodes: return
	while node in edges and len(edges[node]) == 1:
		edge = getone(edges[node])
		assert revnode(edge) in edges
		if len(edges[revnode(edge)]) != 1: break
		if edge[1:] in safe_nodes: break
		safe_nodes.add(edge[1:])
		safe_edges.add(canon(node, edge))
		node = edge

def get_safe_unitigs_and_edges(node_seqs, edges, min_coverage, node_coverage, edge_coverage):
	safe_nodes = set()
	safe_edges = set()
	for node in node_seqs:
		if node not in node_coverage or node_coverage[node] < min_coverage: continue
		add_safe_unitig(">" + node, edges, safe_nodes, safe_edges)
		add_safe_unitig("<" + node, edges, safe_nodes, safe_edges)
	return (safe_nodes, safe_edges)

def remove_and_split_low_coverage(node_seqs, edges, initial_paths, paths_crossing, min_coverage, node_coverage):
	edge_coverage = {}
	for path in initial_paths:
		for i in range(1, len(path)):
			c = canon(path[i-1], path[i])
			if c not in edge_coverage: edge_coverage[c] = 0
			edge_coverage[c] += 1
	(safe_nodes, safe_edges) = get_safe_unitigs_and_edges(node_seqs, edges, min_coverage, node_coverage, edge_coverage)
	remove_edges = set()
	for fromnode in iterate_deterministic(edges):
		for tonode in iterate_deterministic(edges[fromnode]):
			if canon(fromnode, tonode) in safe_edges: continue
			cov = 0
			if canon(fromnode, tonode) in edge_coverage: cov = edge_coverage[canon(fromnode, tonode)]
			if cov >= min_coverage: continue
			remove_edges.add(canon(fromnode, tonode))
			sys.stderr.write("remove low coverage edge " + fromnode + " " + tonode + " with coverage " + str(cov) + "\n")
	for edge in remove_edges:
		(fromnode, tonode) = edge
		assert fromnode in edges
		assert tonode in edges[fromnode]
		edges[fromnode].remove(tonode)
		if (revnode(tonode), revnode(fromnode)) == edge: continue
		assert revnode(tonode) in edges
		assert revnode(fromnode) in edges[revnode(tonode)]
		edges[revnode(tonode)].remove(revnode(fromnode))
	removables = set()
	for node in node_seqs:
		if node in safe_nodes: continue
		if node in node_coverage and node_coverage[node] >= min_coverage: continue
		removables.add(node)
	for path in initial_paths:
		remove_this = False
		add_these = []
		last_break = 0
		if path[0][1:] not in node_seqs or path[0][1:] in removables:
			remove_this = True
			last_break = 1
		for j in range(1, len(path)):
			if path[j][1:] not in node_seqs or path[j][1:] in removables:
				remove_this = True
				if j > last_break: add_these.append(path[last_break:j])
				last_break = j+1
				continue
			(fromnode, tonode) = canon(path[j-1], path[j])
			if fromnode not in edges:
				remove_this = True
				if j > last_break: add_these.append(path[last_break:j])
				last_break = j
				continue
			if tonode not in edges[fromnode]:
				remove_this = True
				if j > last_break: add_these.append(path[last_break:j])
				last_break = j
		if remove_this and last_break < len(path):
			add_these.append(path[last_break:])
		if remove_this: remove_path(paths_crossing, path)
		for addable in add_these:
			if len(addable) >= 2: add_path(paths_crossing, addable)
	for node in iterate_deterministic(removables):
		remove_graph_node(node, node_seqs, edges)
		coverage = 0
		if node in node_coverage: coverage = node_coverage[node]
		sys.stderr.write("removed low coverage node " + node + " with coverage " + str(coverage) + "\n")

def read_node_coverages(node_coverage_file):
	node_coverage = {}
	with open(node_coverage_file) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[0] == 'node': continue
			node_coverage[parts[0]] = float(parts[2])
	return node_coverage

def resolve(node_lens, edge_overlaps, node_seqs, edges, paths_crossing, min_edge_support, min_coverage, removable_nodes):
	maybe_resolvable = set(node_seqs.keys())
	nodes_by_len = [(get_unitig_len(node_lens, edge_overlaps, node_seqs, n), n) for n in maybe_resolvable]
	heapq.heapify(nodes_by_len)
	last_resolved = 0
	while len(nodes_by_len) > 0:
		current_length = nodes_by_len[0][0]
		if current_length > max_resolve_length: break
		current_nodes = set()
		while len(nodes_by_len) > 0 and nodes_by_len[0][0] == current_length:
			(priority, node) = heapq.heappop(nodes_by_len)
			if node not in node_seqs: continue
			current_nodes.add(node)
			nodelen = get_unitig_len(node_lens, edge_overlaps, node_seqs, node)
			assert nodelen == priority
		if len(current_nodes) == 0: continue
		(new_nodes, resolved) = resolve_hairpins(current_length, current_nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, min_edge_support, min_coverage, removable_nodes)
		if len(new_nodes) > 0:
			assert len(resolved) > 0
			sys.stderr.write("resolve k=" + str(current_length) + ", extended " + str(len(resolved)) + " nodes into " + str(len(new_nodes)) + " nodes" + "\n")
			assert current_length > last_resolved
			last_resolved = current_length
			for n in iterate_deterministic(new_nodes):
				if n not in node_seqs: continue # already unitigified
				new_unitig = unitigify_one(node_seqs, node_lens, edges, paths_crossing, n)
				heapq.heappush(nodes_by_len, (get_unitig_len(node_lens, edge_overlaps, node_seqs, new_unitig), new_unitig))
				maybe_resolvable.add(new_unitig)
			continue;
		else:
			assert len(resolved) == 0
		(new_nodes, resolved) = resolve_nodes(current_length, current_nodes, paths_crossing, node_seqs, node_lens, edges, maybe_resolvable, min_edge_support, min_coverage, removable_nodes)
		if len(new_nodes) > 0:
			assert len(resolved) > 0
			sys.stderr.write("resolve k=" + str(current_length) + ", extended " + str(len(resolved)) + " nodes into " + str(len(new_nodes)) + " nodes" + "\n")
		else:
			assert len(resolved) == 0
		if len(resolved) > 0:
			assert current_length > last_resolved
			last_resolved = current_length
		for n in iterate_deterministic(new_nodes):
			if n not in node_seqs: continue # already unitigified
			new_unitig = unitigify_one(node_seqs, node_lens, edges, paths_crossing, n)
			heapq.heappush(nodes_by_len, (get_unitig_len(node_lens, edge_overlaps, node_seqs, new_unitig), new_unitig))
			maybe_resolvable.add(new_unitig)

def get_unitig_len(node_lens, edge_overlaps, node_seqs, node):
	assert node in node_seqs
	assert len(node_seqs[node]) >= 1
	return get_unitig_len_path(node_lens, edge_overlaps, node_seqs[node][0], node_seqs[node][1], node_seqs[node][2])

def get_unitig_len_path(node_lens, edge_overlaps, path, left_clip, right_clip):
	assert len(path) >= 1
	unitig_len = node_lens[path[0][1:]]
	assert node_lens[path[0][1:]] > left_clip
	assert node_lens[path[-1][1:]] > right_clip
	for i in range(1, len(path)):
		overlap = 0
		fromnode = path[i-1]
		tonode = path[i]
		key = canon(fromnode, tonode)
		if key not in edge_overlaps:
			print(path)
			print(key)
		assert key in edge_overlaps
		overlap = edge_overlaps[key]
		if not overlap < node_lens[path[i][1:]]:
			print(overlap)
			print(i)
			print(path)
			print(path[i-1])
			print(path[i])
			print(node_lens[path[i][1:]])
		assert overlap < node_lens[path[i][1:]]
		unitig_len += node_lens[path[i][1:]] - overlap
	assert unitig_len > left_clip + right_clip
	unitig_len -= left_clip + right_clip
	return unitig_len

def iterate_paths(paths_crossing, node):
	paths = [paths_crossing[node][pathid] for pathid in paths_crossing[node]]
	paths.sort()
	for path in paths:
		yield path

def remove_path(paths_crossing, path):
	nodes_in_path = set(n[1:] for n in path)
	for node in nodes_in_path:
		assert node in paths_crossing
		assert id(path) in paths_crossing[node]
		assert paths_crossing[node][id(path)] == path
		assert paths_crossing[node][id(path)] is path
		del paths_crossing[node][id(path)]

def add_path(paths_crossing, path):
	nodes_in_path = set(n[1:] for n in path)
	for node in nodes_in_path:
		if node not in paths_crossing: paths_crossing[node] = {}
		assert id(path) not in paths_crossing[node]
		paths_crossing[node][id(path)] = path

def replace_unitig(node_seqs, node_lens, edges, paths_crossing, unitig):
	global next_unitig_num
	new_node = "unitig_" + str(next_unitig_num)
	next_unitig_num += 1
	node_seqs[new_node] = []
	assert len(unitig) >= 1
	last_node_seq = node_seqs[unitig[0][1:]][0]
	start_clip = node_seqs[unitig[0][1:]][1]
	last_end_clip = node_seqs[unitig[0][1:]][2]
	if unitig[0][0] == "<":
		last_node_seq = [revnode(n) for n in last_node_seq[::-1]]
		start_clip = node_seqs[unitig[0][1:]][2]
		last_end_clip = node_seqs[unitig[0][1:]][1]
	result_node_seq = list(last_node_seq)
	for i in range(1, len(unitig)):
		new_node_seq = node_seqs[unitig[i][1:]][0]
		new_start_clip = node_seqs[unitig[i][1:]][1]
		new_end_clip = node_seqs[unitig[i][1:]][2]
		if unitig[i][0] == "<":
			new_node_seq = [revnode(n) for n in new_node_seq[::-1]]
			new_start_clip = node_seqs[unitig[i][1:]][2]
			new_end_clip = node_seqs[unitig[i][1:]][1]
		overlap = edge_overlaps[canon(unitig[i-1], unitig[i])]
		found = False
		for j in range(0, len(last_node_seq)):
			if new_start_clip >= node_lens[last_node_seq[j][1:]]: continue
			if get_unitig_len_path(node_lens, edge_overlaps, last_node_seq[j:], new_start_clip, last_end_clip) == overlap:
				assert last_node_seq[j:] == new_node_seq[:len(last_node_seq)-j]
				result_node_seq += new_node_seq[len(last_node_seq)-j:]
				found = True
				break
		if not found:
			key = canon(last_node_seq[-1], new_node_seq[0])
			if key not in edge_overlaps:
				print(unitig)
				print(key)
				print(last_end_clip)
				print(new_start_clip)
				print(overlap)
				print(last_node_seq)
				print(new_node_seq)
				print(node_seqs[unitig[i-1][1:]])
				print(node_seqs[unitig[i][1:]])
			assert key in edge_overlaps
			if not edge_overlaps[key] == overlap + last_end_clip + new_start_clip:
				print(unitig)
				print(key)
				print(last_end_clip)
				print(new_start_clip)
				print(overlap)
				print(last_node_seq)
				print(new_node_seq)
				print(node_seqs[unitig[i-1][1:]])
				print(node_seqs[unitig[i][1:]])
			assert edge_overlaps[key] == overlap + last_end_clip + new_start_clip
			result_node_seq += new_node_seq
		last_node_seq = new_node_seq
		last_end_clip = new_end_clip
	node_seqs[new_node] = (result_node_seq, start_clip, last_end_clip)
	for i in range(1, len(node_seqs[new_node][0])):
		key = canon(node_seqs[new_node][0][i-1], node_seqs[new_node][0][i])
		if key not in edge_overlaps:
			print(unitig)
			for node in unitig: print(node_seqs[node[1:]])
			print(node_seqs[new_node])
			print(key)
		assert key in edge_overlaps
	edges[">" + new_node] = set()
	edges["<" + new_node] = set()
	if unitig[-1] in edges:
		add_edges = []
		for edge in iterate_deterministic(edges[unitig[-1]]):
			assert revnode(edge) in edges
			add_edges.append((">" + new_node, edge))
			edge_overlaps[canon(">" + new_node, edge)] = edge_overlaps[canon(unitig[-1], edge)]
		for edge in add_edges:
			edges[edge[0]].add(edge[1])
			edges[revnode(edge[1])].add(revnode(edge[0]))
	if revnode(unitig[0]) in edges:
		add_edges = []
		for edge in iterate_deterministic(edges[revnode(unitig[0])]):
			assert revnode(edge) in edges
			add_edges.append(("<" + new_node, edge))
			edge_overlaps[canon("<" + new_node, edge)] = edge_overlaps[canon(revnode(unitig[0]), edge)]
		for edge in add_edges:
			edges[edge[0]].add(edge[1])
			edges[revnode(edge[1])].add(revnode(edge[0]))
	if not len(unitig) == len(set(n[1:] for n in unitig)):
		print(unitig)
		for n in unitig: print(node_seqs[n[1:]])
	assert len(unitig) == len(set(n[1:] for n in unitig))
	for node in unitig:
		remove_graph_node(node[1:], node_seqs, edges)
	replaceable_paths = []
	found_identities = set()
	for node in unitig:
		for path in iterate_paths(paths_crossing, node[1:]):
			if id(path) in found_identities: continue
			found_identities.add(id(path))
			replaceable_paths.append(path)
	for path in replaceable_paths:
		remove_path(paths_crossing, path)
	is_forward = {}
	for n in unitig:
		if n[0] == ">":
			is_forward[n[1:]] = True
		else:
			is_forward[n[1:]] = False
	unitig_nodes = set(n[1:] for n in unitig)
	assert len(unitig_nodes) == len(unitig)
	for path in replaceable_paths:
		new_path = []
		next_start = 0
		for i in range(0, len(path)):
			if i < next_start: continue
			if path[i][1:] not in unitig_nodes:
				new_path.append(path[i])
				continue
			path_end = i+1
			while path_end < len(path) and path[path_end][1:] in unitig_nodes and not (path[path_end] == unitig[0] and path[path_end-1] == unitig[-1]) and not (path[path_end] == revnode(unitig[-1]) and path[path_end-1] == revnode(unitig[0])):
				path_end += 1
			next_start = path_end
			new_path.append((">" if (is_forward[path[i][1:]] == (path[i][0] == ">")) else "<") + new_node)
		add_path(paths_crossing, new_path)
	return new_node

def extend_forward(node, edges):
	result = [node]
	while True:
		pos = result[-1]
		if pos not in edges: break
		if len(edges[pos]) != 1: break
		if len(edges[pos]) == 1 and getone(edges[pos])[1:] == pos[1:]: break # palindrome hairpin
		newpos = getone(edges[pos])
		assert revnode(newpos) in edges
		if len(edges[revnode(newpos)]) != 1: break
		if len(edges[pos]) == 1 and getone(edges[pos]) == node: # circular unitig
			result.append(node)
			break
		if len(edges[pos]) == 1 and getone(edges[pos])[1:] == node[1:]: break # circular unitig
		result.append(newpos)
	return result

def unitigify_one(node_seqs, node_lens, edges, paths_crossing, node):
	forward_extension = extend_forward(">" + node, edges)
	assert len(forward_extension) >= 1
	backward_extension = []
	if len(forward_extension) == 1 or forward_extension[0] != forward_extension[-1]:
		backward_extension = extend_forward("<" + node, edges)
		assert len(backward_extension) >= 1
	if len(forward_extension) + len(backward_extension) == 2: return node
	unitig = [revnode(n) for n in backward_extension[::-1]][:-1] + forward_extension
	if unitig[0] == unitig[-1]: unitig.pop()
	# sys.stderr.write("unitigify " + str(unitig) + " (start " + node + ", " + str(forward_extension) + "," + str(backward_extension) + ")")
	new_name = replace_unitig(node_seqs, node_lens, edges, paths_crossing, unitig)
	# sys.stderr.write(" into " + new_name + "\n")
	return new_name

def unitigify_all(node_seqs, node_lens, edges, paths_crossing):
	maybe_unitigifiable = set(node_seqs)
	for node in iterate_deterministic(maybe_unitigifiable):
		if node not in node_seqs: continue # already unitigified
		unitigify_one(node_seqs, node_lens, edges, paths_crossing, node)

def get_seq(base_seqs, edge_overlaps, nodeseq, left_clip, right_clip):
	assert len(nodeseq) >= 1
	assert left_clip >= 0
	assert right_clip >= 0
	result = base_seqs[nodeseq[0][1:]]
	if nodeseq[0][0] == "<": result = revcomp(result)
	for i in range(1, len(nodeseq)):
		add_seq = base_seqs[nodeseq[i][1:]]
		if nodeseq[i][0] == "<": add_seq = revcomp(add_seq)
		key = canon(nodeseq[i-1], nodeseq[i])
		assert key in edge_overlaps
		overlap = edge_overlaps[key]
		result += add_seq[overlap:]
	assert len(result) > left_clip + right_clip
	if left_clip > 0: result = result[left_clip:]
	if right_clip > 0: result = result[:-right_clip]
	return result

(base_seqs, node_seqs, edges, edge_overlaps) = read_graph(input_gfa)

node_lens = {}
paths_crossing = {}
for node in base_seqs:
	node_lens[node] = len(base_seqs[node])
	paths_crossing[node] = {}

del base_seqs

initial_paths = []

for l in sys.stdin:
	parts = l.strip().split('\t')
	if len(parts) == 1:
		path = parts[0].replace(">", "\t>").replace("<", "\t<").strip().split('\t')
	else:
		pathstr = parts[5]
		path = pathstr.replace(">", "\t>").replace("<", "\t<").strip().split('\t')
	if len(path) <= 1: continue
	add_path(paths_crossing, path)
	initial_paths.append(path)

node_coverage = read_node_coverages(node_coverage_file)
removable_nodes = set()
for n in node_lens:
	if n not in node_coverage or node_coverage[n] < min_allowed_coverage:
		removable_nodes.add(n)
remove_and_split_low_coverage(node_seqs, edges, initial_paths, paths_crossing, 0, node_coverage)
del initial_paths
unitigify_all(node_seqs, node_lens, edges, paths_crossing)

for coverage in resolve_steps:
	sys.stderr.write("resolve with edge support " + str(coverage) + "\n")
	resolve(node_lens, edge_overlaps, node_seqs, edges, paths_crossing, coverage, min_allowed_coverage, removable_nodes)

del node_lens

sys.stderr.write("done resolving" + "\n")

unitig_name = {}
unitig_num = 1
for n in iterate_deterministic(node_seqs):
	if len(node_seqs[n]) == 1:
		assert node_seqs[n][0][0] == ">"
		unitig_name[n] = node_seqs[n][0][1:]
	else:
		unitig_name[n] = "unitig_" + str(unitig_num) + "_" + "_".join(n[1:] + "n" + ("f" if n[0] == ">" else "b") for n in node_seqs[n][0])
		unitig_num += 1

sys.stderr.write("write paths" + "\n")

found_identities = set()
with open(out_path_file, "w") as f:
	for n in iterate_deterministic(paths_crossing):
		for path in iterate_paths(paths_crossing, n):
			if id(path) in found_identities: continue
			found_identities.add(id(path))
			f.write("".join(n[0] + unitig_name[n[1:]] for n in path) + "\n")

del paths_crossing
del found_identities

sys.stderr.write("paths written" + "\n")
sys.stderr.write("write graph" + "\n")

for e1 in iterate_deterministic(edges):
	for e2 in iterate_deterministic(edges[e1]):
		print("L\t" + unitig_name[e1[1:]] + "\t" + ("+" if e1[0] == ">" else "-") + "\t" + unitig_name[e2[1:]] + "\t" + ("+" if e2[0] == ">" else "-") + "\t" + str(edge_overlaps[canon(e1, e2)]) + "M")

del edges

base_seqs = read_graph_only_bases(input_gfa)

for n in iterate_deterministic(node_seqs):
	print("S\t" + unitig_name[n] + "\t" + get_seq(base_seqs, edge_overlaps, node_seqs[n][0], node_seqs[n][1], node_seqs[n][2]))

sys.stderr.write("graph written" + "\n")
sys.stderr.write("write name mapping" + "\n")

with open(resolve_namemapping_file, "w") as f:
	for n in node_seqs:
		f.write(unitig_name[n] + "\t" + "".join(node_seqs[n][0]) + ":" + str(node_seqs[n][1]) + ":" + str(node_seqs[n][2]) + "\n")

sys.stderr.write("name mapping written" + "\n")
