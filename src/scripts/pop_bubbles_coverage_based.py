#!/usr/bin/env python

import sys
import graph_functions as gf
import networkx as nx

graph_file = sys.argv[1]
node_coverage_file = sys.argv[2]
haploid = gf.str2bool(sys.argv[3])

# gfa to stdout

max_bubble_pop_size = 10
max_poppable_node_size = 200000
min_coverage_node_size = max_poppable_node_size / 2
max_poppable_coverage = 0
max_coverage_delta = 1.5

if haploid: max_coverage_delta += 0.3

def find_component_coverage(key, belongs_to_component, component_coverage_sum, component_length_sum):
	comp_coverage = avg_coverage
	if key in belongs_to_component:
		if belongs_to_component[key] in component_coverage_sum:
			if component_long_nodes[belongs_to_component[key]] > 5:
				comp_coverage = component_coverage_sum[belongs_to_component[key]] / component_length_sum[belongs_to_component[key]]
	return comp_coverage

def merge(parent, chain_coverage_sum, chain_length_sum, left, right):
	left = gf.find(parent, left)
	right = gf.find(parent, right)
	assert parent[left] == left
	assert parent[right] == right
	assert parent[left] != parent[right]
	chain_coverage_sum[left] += chain_coverage_sum[right]
	chain_length_sum[left] += chain_length_sum[right]
	parent[right] = left

def remove_graph_node(node, edges):
	if ">" + node in edges:
		for edge in edges[">" + node]:
			assert gf.revnode(edge) in edges
			assert "<" + node in edges[gf.revnode(edge)]
			edges[gf.revnode(edge)].remove("<" + node)
		del edges[">" + node]
	if "<" + node in edges:
		for edge in edges["<" + node]:
			assert gf.revnode(edge) in edges
			assert ">" + node in edges[gf.revnode(edge)]
			edges[gf.revnode(edge)].remove(">" + node)
		del edges["<" + node]

# Detecting Superbubbles in Assembly Graphs, Onodera et al 2013
# fig. 5
def find_bubble(s, edges, nodelens, max_size):
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
		if v != s and nodelens[v[1:]] > max_size: return None
		if v not in edges: return None
		if len(edges[v]) == 0: return None
		for u in edges[v]:
			if u[1:] == v[1:]: return None
			if gf.revnode(u) in visited: return None
			if u == s: return None
			assert u not in visited
			seen.add(u)
			assert gf.revnode(u) in edges
			assert len(edges[gf.revnode(u)]) >= 1
			has_nonvisited_parent = False
			for parent_edge in edges[gf.revnode(u)]:
				parent = gf.revnode(parent_edge)
				if parent not in visited: has_nonvisited_parent = True
			if not has_nonvisited_parent: S.append(u)
		if len(S) == 1 and len(seen) == 1 and S[0] == gf.getone(seen):
			t = S.pop()
			if t in edges:
				for edge in edges[t]:
					if edge == s: return None
			return (s, t)
	return None

def pop_bubble(start, end, removed_nodes, removed_edges, edges, coverage, nodelens, conservative = False):
	bubble_nodes = set()
	bubble_edges = set()
	visited	  = set()
	max_bubble_node_size = 0
	predecessor = {}
	stack = []
	stack.append((start, start, coverage.get(start[1:], 0)))
	visited.add(start)

	while len(stack) > 0:
		(top, before, pathwidth) = stack.pop()
		if top in predecessor:
			p = predecessor[top]
		else:
			p = ""
		if top not in predecessor: predecessor[top] = (before, pathwidth)
		# special case for three-node bubbles (e.g. transitive), we won't replace the predecessor to skip the intermediate node unless it has really low coverage
		# without this we'd force keep the begin node to end node edge if the end node has higher coverage than the intermediate node, even if they are 22.1 vs 22.09
		if (top != end or before != start or predecessor[top][1] < 0.5*avg_coverage) and predecessor[top][1] < pathwidth: predecessor[top] = (before, pathwidth)
		bubble_nodes.add(top[1:])
		bubble_edges.add((before, top))
		if top[1:] != start[1:] and top[1:] != end[1:] and nodelens[top[1:]] > max_bubble_node_size:
			max_bubble_node_size = nodelens[top[1:]]
		if top == end: continue
		for edge in gf.iterate_deterministic(edges[top], end):
			if edge not in visited:
				visited.add(edge)
				stack.append((edge, top, min(pathwidth, coverage.get(edge[1:], 0))))
	assert end in predecessor
	#sys.stderr.write("Processing bubble from %s to %s and conservative mode is %s\n"%(start, end, conservative))
	# TODO: when our coverage is more trustworthy (e.g. low coverage isn't due to no unique nodes in a path, we can run this in conservative popping mode to remove noise, until then do nothing in these bubbles, also do nothign for 3-node bubbles in high coverage
	if (conservative == True and (len(bubble_nodes) == 3 or len(bubble_nodes) > max_bubble_pop_size / 2)) or len(bubble_nodes) > max_bubble_pop_size:
		#sys.stderr.write("Error bubble beteween %s and %s is not poppable because it its size of %s is  larger than max %s\n"%(start, end, len(bubble_nodes), max_bubble_pop_size))
		return

	path = [end]
	while path[-1] != start:
		path.append(predecessor[path[-1]][0])
	path = path[::-1]
	kept_nodes = set()
	kept_edges = set()
	for node in path:
		kept_nodes.add(node[1:])
	for i in range(1, len(path)):
		kept_edges.add((path[i-1], path[i]))
	assert len(kept_nodes) == len(path)
	assert len(kept_edges) == len(path)-1
	assert len(kept_nodes) <= len(bubble_nodes)
	assert len(kept_edges) < len(bubble_edges) or len(bubble_edges) == 1

	# check that the length is acceptable for standard bubbles, transitive edge bubbles are allowed a larger pop length
	if len(bubble_nodes) > 3 and max_bubble_node_size > max_poppable_node_size: return

	# set minimum coverage, for 3-node (transitive) bubbles we are conservative always, otherwise be agressive in haploid genomes and when we are surrounded by large likely resolved nodes (e.g. within haplotype bubble)
	comp_coverage = find_component_coverage(start, belongs_to_component, component_coverage_sum, component_length_sum)
	if conservative == True:
		max_poppable_coverage = 0.25 * comp_coverage
	elif start[1:] not in coverage or end[1:] not in coverage:
		max_poppable_coverage = 0.25 * comp_coverage
	elif len(bubble_nodes) == 3:
		max_poppable_coverage = 0.5*comp_coverage
	elif start[1:] in coverage and end[1:] in coverage and nodelens[start[1:]] > max_poppable_node_size and nodelens[end[1:]] > max_poppable_node_size and coverage[start[1:]] <= 1.5*comp_coverage and coverage[start[1:]] >= 0.5 * comp_coverage and coverage[end[1:]] <= 1.5*comp_coverage and coverage[end[1:]] >= 0.5*comp_coverage:
		max_poppable_coverage = comp_coverage
	elif haploid:
		max_poppable_coverage = comp_coverage
	else:
		max_poppable_coverage = 0.5*comp_coverage
	#sys.stderr.write("Processing bubble from %s to %s and nodes %s and max poppable coverage is %s\n"%(start, end, bubble_nodes, max_poppable_coverage))

	for node in bubble_nodes:
		if node in kept_nodes: continue
		c = (coverage[node] if node in coverage else 0)
		if c > max_poppable_coverage:
			kept_nodes.add(node)
			continue
		remove_graph_node(node, edges)
		removed_nodes.add(node)
	for edge in bubble_edges:
		if edge in kept_edges: continue
		if (gf.revnode(edge[1]), gf.revnode(edge[0])) in kept_edges: continue
		# if we have a bubble where there is a connection between the start and end and we decided to keep a node in between them, remove the edge skipping that node then
		if set([start,end]) == kept_nodes or edge[0] != start or edge[1] != end:
			if edge[0][1:] in kept_nodes or edge[1][1:] in kept_nodes: continue
		removed_edges.add(edge)
		if edge[0] in edges:
			if edge[1] in edges[edge[0]]:
				edges[edge[0]].remove(edge[1])
		if gf.revnode(edge[1]) in edges:
			if gf.revnode(edge[0]) in edges[gf.revnode(edge[1])]:
				edges[gf.revnode(edge[1])].remove(gf.revnode(edge[0]))

def try_pop_tip(start, edges, coverage, removed_nodes, removed_edges, max_removable, nodelens):
	if start not in edges: return
	#sys.stderr.write("Processing pop tip from node %s with edges %s max removable is %s\n"%(start, edges[start], max_removable))
	if len(edges[start]) < 2 or len(edges[start]) > 4: return
	max_coverage = 0
	max_len = 0
	keeps = []
	for node in gf.iterate_deterministic(edges[start]):
		assert gf.revnode(node) in edges
		if len(edges[gf.revnode(node)]) != 1: keeps.append(node)
		if node in edges and len(edges[node]) > 0: keeps.append(node)
		coverage_here = 0
		if node[1:] in coverage: coverage_here = coverage[node[1:]]
		if coverage_here > max_coverage:
			max_coverage = coverage_here
		if nodelens[node[1:]] > max_len:
			max_len = nodelens[node[1:]]
	remove_this = []
	for node in gf.iterate_deterministic(edges[start]):
		if node in keeps:
			continue
		coverage_here = 0
		if node[1:] in coverage: coverage_here = coverage[node[1:]]
		# if we have a node w/better coverage, remove this one. Ties are broken by length
		if coverage_here < max_removable and len(remove_this)+1 < len(edges[start]) and (coverage_here < max_coverage or (coverage_here == max_coverage and nodelens[node[1:]] < max_len)):
			remove_this.append(node)
	if len(remove_this) == 0: return
	for remove_node in remove_this:
		assert remove_node[1:] in nodelens
		if nodelens[remove_node[1:]] > min_coverage_node_size: continue
		removed_nodes.add(remove_node[1:])
		removed_edges.add((start, remove_node))

coverage = {}
with open(node_coverage_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "node" and parts[1] == "coverage": continue
		coverage[parts[0]] = float(parts[1])

nodelens = {}
edges = {}
nodelines = []
edgelines = []
removed_nodes = set()
removed_edges = set()

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodelines.append((parts[1], l.strip()))
			nodelens[parts[1]] = len(parts[2])
		elif parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			edgelines.append((fromnode, tonode, l.strip()))
			if fromnode not in edges: edges[fromnode] = set()
			if gf.revnode(tonode) not in edges: edges[gf.revnode(tonode)] = set()
			edges[fromnode].add(tonode)
			edges[gf.revnode(tonode)].add(gf.revnode(fromnode))

long_coverage_cov_sum = 0.0
long_coverage_len_sum = 0.0
for node in nodelens:
	if node not in coverage: continue
	if nodelens[node] < min_coverage_node_size: continue
	long_coverage_len_sum += nodelens[node]
	long_coverage_cov_sum += nodelens[node] * coverage[node]

avg_coverage = 0
if long_coverage_len_sum != 0:
	avg_coverage = long_coverage_cov_sum / long_coverage_len_sum
component_coverage_sum = {}
component_length_sum = {}
component_long_nodes = {}
belongs_to_component = {}
# we load the graph twice rather than re-implementing connected component search
G = nx.Graph()
gf.load_indirect_graph(graph_file, G)
comp = 0
for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
	comp += 1

	for node in current_component:
		belongs_to_component[node] = comp
		if nodelens[node] < min_coverage_node_size: continue

		if comp not in component_coverage_sum:
			component_coverage_sum[comp] = 0
			component_length_sum[comp] = 0
			component_long_nodes[comp] = 0
		if node in coverage:
			component_length_sum[comp] += nodelens[node]
			component_coverage_sum[comp] += coverage[node] * nodelens[node]
			component_long_nodes[comp] += 1

chain_coverage_sum = {}
chain_length_sum = {}
parent = {}
for node in nodelens:
	parent[node] = node
	if node in coverage:
		chain_length_sum[node] = nodelens[node]
		chain_coverage_sum[node] = coverage[node] * nodelens[node]
	else:
		chain_length_sum[node] = nodelens[node]
		chain_coverage_sum[node] = 0

possible_merges = []

for edge in edges:
	bubble = find_bubble(edge, edges, nodelens, max_poppable_node_size * 5)
	if not bubble: continue
	cov = max(coverage[bubble[0][1:]] / coverage[bubble[1][1:]], coverage[bubble[1][1:]] / coverage[bubble[0][1:]]) if bubble[0][1:] in coverage and bubble[1][1:] in coverage else 0
	possible_merges.append((bubble[0][1:], bubble[1][1:], cov))

possible_merges.sort(key=lambda x: x[2])
while True:
	new_possible_merges = []
	merged_any = False
	for triplet in possible_merges:
		(node1, node2, difference) = triplet
		key1 = gf.find(parent, node1)
		key2 = gf.find(parent, node2)
		if key1 == key2: continue
		chain1_cov = chain_coverage_sum[gf.find(parent, node1)] / chain_length_sum[gf.find(parent, node1)]
		chain2_cov = chain_coverage_sum[gf.find(parent, node2)] / chain_length_sum[gf.find(parent, node2)]
		assert(belongs_to_component[key1] == belongs_to_component[key2])
		comp_coverage = find_component_coverage(key1, belongs_to_component, component_coverage_sum, component_length_sum)
		if chain1_cov > chain2_cov * max_coverage_delta and (chain1_cov > comp_coverage * 1.5 or chain2_cov > comp_coverage * 1.5):
			new_possible_merges.append(triplet)
			continue
		if chain2_cov > chain1_cov * max_coverage_delta and (chain1_cov > comp_coverage * 1.5 or chain2_cov > comp_coverage * 1.5):
			new_possible_merges.append(triplet)
			continue
		merge(parent, chain_coverage_sum, chain_length_sum, node1, node2)
		merged_any = True
	if not merged_any: break
	possible_merges = new_possible_merges

unique_chains = set()
tip_chains = set()
for node in nodelens:
	key = gf.find(parent, node)
	if key not in chain_coverage_sum: continue
	comp_coverage = find_component_coverage(key, belongs_to_component, component_coverage_sum, component_length_sum)
	chain_coverage = chain_coverage_sum[key] / chain_length_sum[key]
	#sys.stderr.write("Checking chain with key %s and coverqge %s versus avg coverage %s and componenet %s\n"%(key, chain_coverage, avg_coverage, comp_coverage))
	if chain_coverage <= comp_coverage * 2.5:
		tip_chains.add(key)
	if chain_coverage >= comp_coverage * 0.5 and chain_coverage <= comp_coverage * 2.5:
		unique_chains.add(key)

for node in gf.iterate_deterministic(nodelens):
	key = gf.find(parent, node)
	if key not in unique_chains: continue
	if node in removed_nodes: continue
	comp_coverage = find_component_coverage(key, belongs_to_component, component_coverage_sum, component_length_sum)
	chain_coverage = chain_coverage_sum[key] / chain_length_sum[key]

	bubble = find_bubble(">" + node, edges, nodelens, max_poppable_node_size * 5)
	if bubble:
		assert bubble[0] == ">" + node
		assert bubble[1][1:] != node
		if gf.find(parent, bubble[1][1:]) != key: continue
		pop_bubble(bubble[0], bubble[1], removed_nodes, removed_edges, edges, coverage, nodelens, chain_coverage > comp_coverage * 1.5)
	bubble = find_bubble("<" + node, edges, nodelens, max_poppable_node_size * 5)
	if bubble:
		assert bubble[0] == "<" + node
		assert bubble[1][1:] != node
		if gf.find(parent, bubble[1][1:]) != key: continue
		pop_bubble(bubble[0], bubble[1], removed_nodes, removed_edges, edges, coverage, nodelens, chain_coverage > comp_coverage * 1.5)

for node in gf.iterate_deterministic(nodelens):
	key = gf.find(parent, node)
	if key not in tip_chains: continue
	if node in removed_nodes: continue
	bubble = find_bubble(">" + node, edges, nodelens, max_poppable_node_size)
	comp_coverage = find_component_coverage(node, belongs_to_component, component_coverage_sum, component_length_sum)
	chain_coverage = chain_coverage_sum[key] / chain_length_sum[key]

	if not bubble:
		#sys.stderr.write("Pop tip called with chain cov %s and comp %s\n"%(chain_coverage, comp_coverage))
		try_pop_tip(">" + node, edges, coverage, removed_nodes, removed_edges, comp_coverage * 0.75 if chain_coverage <= comp_coverage * 1.5 else comp_coverage * 0.25, nodelens)
	bubble = find_bubble("<" + node, edges, nodelens, max_poppable_node_size)
	if not bubble:
		#sys.stderr.write("Pop tip reverse called with chain cov %s and comp%s\n"%(chain_coverage, comp_coverage))
		try_pop_tip("<" + node, edges, coverage, removed_nodes, removed_edges, comp_coverage * 0.75 if chain_coverage <= comp_coverage * 1.5 else comp_coverage * 0.25, nodelens)

for node in removed_nodes:
	sys.stderr.write(node + "\n")

for edge in removed_edges:
	if edge[0][1:] in removed_nodes: continue
	if edge[1][1:] in removed_nodes: continue
	sys.stderr.write(edge[0] + "\t" + edge[1] + "\n")

for node in nodelines:
	if node[0] in removed_nodes: continue
	print(node[1])

for edge in edgelines:
	if edge[0][1:] in removed_nodes: continue
	if edge[1][1:] in removed_nodes: continue
	if (edge[0], edge[1]) in removed_edges: continue
	if (gf.revnode(edge[1]), gf.revnode(edge[0])) in removed_edges: continue
	print(edge[2])
