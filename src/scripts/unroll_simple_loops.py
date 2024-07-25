#!/usr/bin/env python

import sys
import graph_functions as gf

min_len = 100000
graph_file = sys.argv[1]
cov_file = sys.argv[2]
# graph to stdout

nodelens = {}
coverage = {}
long_coverage_len_sum = 0.0
long_coverage_cov_sum = 0.0
edge_overlaps = {}
edges = {}

with open(cov_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if "node" in parts[0]: continue
		coverage[parts[0]] = float(parts[1])
		nodelens[parts[0]] = int(parts[2])
		if nodelens[parts[0]] > min_len:
			long_coverage_len_sum += nodelens[parts[0]]
			long_coverage_cov_sum += nodelens[parts[0]] * coverage[parts[0]]

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodelens[parts[1]] = len(parts[2])
		if parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromnode not in edge_overlaps:
				edge_overlaps[fromnode] = set()
			edge_overlaps[fromnode].add(tonode)
			edges[gf.canontip(fromnode, tonode)] = int(parts[5][:-1])

avg_coverage = long_coverage_cov_sum / long_coverage_len_sum

tounroll = set()
copy2edges = set()
removed = set()
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if ">" + parts[1] not in edge_overlaps or "<" + parts[1] not in edge_overlaps: continue
			if len(edge_overlaps[">" + parts[1]]) > 1 or len(edge_overlaps["<" + parts[1]]) > 1: continue
			tonodeFwd = next(iter(edge_overlaps[">" + parts[1]]))
			tonodeRev = next(iter(edge_overlaps["<" + parts[1]]))
			#sys.stderr.write("Found neighbors %s and %s and selected %s and %s with set %s\n"%(edge_overlaps[">" + parts[1]], edge_overlaps[">" + parts[1]], tonodeFwd, tonodeRev, set([tonodeFwd[:1], tonodeRev[:1]])))
			if tonodeFwd in edge_overlaps[tonodeRev] or tonodeRev in edge_overlaps[tonodeFwd] or tonodeFwd[1:] == parts[1] or tonodeFwd[1:] != tonodeRev[1:] or len(set([tonodeFwd[:1], tonodeRev[:1]])) < 2: continue
			#sys.stderr.write("Checking loop at %s with node %s with coverage of %s and %s and theshold 0.5 is %s and 2.5 is %s\n"%(tonodeFwd, parts[1], (coverage[parts[1]] if parts[1] in coverage else "NA"), (coverage[tonodeFwd[1:]] if tonodeFwd[1:] in coverage else"NA"), 0.5*avg_coverage, 2.5*avg_coverage)) 
			if parts[1] not in coverage or tonodeFwd[1:] not in coverage or coverage[parts[1]] > 1.5 * avg_coverage or coverage[tonodeFwd[1:]] <= 0.5 * avg_coverage or coverage[tonodeFwd[1:]] >= 2.5 * avg_coverage:
				# if the loop has no coverage and is contained within its overlap and the doubly traversed node has single copy coverage has single copy count, we remove the loop by removing a node not unrolling
				if parts[1] not in coverage and tonodeFwd[1:] in coverage and coverage[tonodeFwd[1:]] <= 1.5 * avg_coverage and nodelens[parts[1]] - edges[gf.canontip(">"+parts[1], tonodeFwd)] < 10:
					removed.add(parts[1])
					continue
				else:
 					continue
			# if we already unrolled something from this node don't do it again, we can do it in the next round
			if tonodeFwd[1:] in tounroll:
				continue
			neighbors = edge_overlaps[">" + tonodeFwd[1:]].union(edge_overlaps["<" + tonodeFwd[1:]])
			neighbors.remove(">" + parts[1])
			neighbors.remove("<" + parts[1])
			if len(neighbors) == 0: continue
			skip=False
			# check if the neighbors are acceptable, if the coverage is too high (we use 2 in case of nested loop) or if the nodes are very short, don't resolve
			if len(neighbors) <= 2:
				for n in neighbors:
					#sys.stderr.write("Checking neighbors for node %s which is %s and it has coverage %s and length %s vs min %s\n"%(parts[1], n[1:], coverage[n[1:]], nodelens[n[1:]], min_len))
					# skip the unroll when we have too high of a coverage or the length of the surrounding nodes is too short or we have a self loop (n[1:] == tonodeFwd) or we already unrolled a neighbor (n[1:] in tounroll)
					if (n[1:] in coverage and coverage[n[1:]] >= 2.5 * avg_coverage) or nodelens[n[1:]] < int(min_len / 6) or n[1:] == tonodeFwd[1:] or n[1:] in tounroll:
						skip=True
			else:
				skip=True
			if skip: continue
			#sys.stderr.write("Found candidate loop with node %s and the sink is %s and %s and neighbors are %s\n"%(parts[1], tonodeFwd, tonodeRev, neighbors))
			neighbor = neighbors.pop()
			# we will duplicate this node
			tounroll.add(tonodeFwd[1:])
			sys.stderr.write("%s%s\t>%s:0:0\n"%(tonodeFwd[1:], "_copy1", tonodeFwd[1:]))
			sys.stderr.write("%s%s\t>%s:0:0\n"%(tonodeFwd[1:], "_copy2", tonodeFwd[1:]))
			# add the two sides of the repeat node to the set so we know which edges to update to point to the newly created duplicate node
			copy2edges.add(neighbor)
			second_copy = ""
			for s in iter(edge_overlaps[neighbor]):
				if s[1:] == tonodeFwd[1:]:
					second_copy = gf.revnode(s)
					break 
			assert(second_copy != "")
			#sys.stderr.write("Neighbors are %s of second copy node %s and will be getting next node %s\n"%(neighbors, second_copy, edge_overlaps[neighbor]))
			copy2edges.add(next(iter(edge_overlaps[second_copy].difference(neighbors))))

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			if parts[1] in tounroll:
				print(parts[0] + "\t" + parts[1] + "_copy1" + "\t" + "\t".join(parts[2:]))
				print(parts[0] + "\t" + parts[1] + "_copy2" + "\t" + "\t".join(parts[2:]))
			elif parts[1] in removed:
				continue
			else:
				print(l.strip())
		if parts[0] == 'L':
			if parts[1] in removed or parts[3] in removed: continue

			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			if (fromnode in copy2edges and parts[3] in tounroll):
				print("\t".join(parts[0:4]) + "_copy2" + "\t" + "\t".join(parts[4:]))
			elif (tonode in copy2edges and parts[1] in tounroll):
				print("\t".join(parts[0:2]) + "_copy2" + "\t" + "\t".join(parts[2:]))
			elif parts[1] in tounroll:
				print("\t".join(parts[0:2]) + "_copy1" + "\t" + "\t".join(parts[2:]))
			elif parts[3] in tounroll:
				print("\t".join(parts[0:4]) + "_copy1" + "\t" + "\t".join(parts[4:]))
			else: 
				print(l.strip())
