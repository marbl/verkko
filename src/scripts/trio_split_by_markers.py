#!/usr/bin/env python

import sys

MIN_DISTANCE=1000000

# this reads a bed file with the 4th column specifying the haplotype, assumes blocks have been reasonably merged to avoid noise
def read_hap(fin, nodes):
	prev=""
	prev_hap = ""
	srt = 2 ** 32 - 1
	end = 0
	with open(fin) as f:
		for l in f:
			parts = l.strip().split('\t')
			if parts[3] not in nodes: nodes[parts[3]] = {}

			# when we hit a new node or a new haplotype, we record the current interval and start a new one
			if prev != parts[0] or prev_hap != parts[3]:
				if prev != parts[0] and prev in nodelens and end < nodelens[prev]: end = nodelens[prev]
				if prev != "" and end - srt >= MIN_DISTANCE:
					#sys.stderr.write("Recording region %s-%s in %s with hap %s\n"%(srt, end, prev, prev_hap))
					if prev not in nodes[prev_hap]: nodes[prev_hap][prev] = []
					nodes[prev_hap][prev].append(srt)
					nodes[prev_hap][prev].append(end)
				end = 0
				srt = 2**32-1
			prev = parts[0]
			prev_hap = parts[3]
			if srt > int(parts[1]): srt = int(parts[1])
			if end < int(parts[2]): end = int(parts[2])
		if prev in nodelens and end < nodelens[prev]: end = nodelens[prev]	# if we don't have any markers between last and end assume it's the same haplotype as we last saw
		if end - srt >= MIN_DISTANCE:
			#sys.stderr.write("Recording region %s-%s in %s with hap %s\n"%(srt, end, prev, prev_hap))
			if prev not in nodes[prev_hap]: nodes[prev_hap][prev] = []
			nodes[prev_hap][prev].append(srt)
			nodes[prev_hap][prev].append(end)

in_graph_file = sys.argv[1]
markers_file = sys.argv[2]
out_mapping_file = sys.argv[3]

nodelens = {}
with open(in_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'S':
			nodelens[parts[1]] = len(parts[2])

haps = {}
read_hap(markers_file, haps)
cut_positions = {}

for node in nodelens:
	# combine breaks from all the haplotypes together for this node
	breaks=[]
	numhaps = 0
	for hap in haps:
		if node not in haps[hap]: continue
		breaks += haps[hap][node]
		numhaps += 1

	# when we only have one hap in a node, nothing to do
	if numhaps <= 1: continue
	sys.stderr.write("Found candidate node %s of len %s with breaks %s\n"%(node, nodelens[node], breaks))

    # sort the positions of switches and try to break the node as needed
	breaks=sorted(breaks)
	last_cut = 0
	for b in breaks:
		if b < MIN_DISTANCE or b + MIN_DISTANCE > nodelens[node] or last_cut + MIN_DISTANCE > b: continue # too close to start/end/other break
		if node not in cut_positions: cut_positions[node] = []
		cut_positions[node].append(b)
		last_cut = b
		sys.stderr.write("Breaking node %s of len %s at %s\n"%(node, nodelens[node], b))

with open(out_mapping_file, "a") as f:
	for node in cut_positions:
		if len(cut_positions[node]) > 0:
			for i in range(0, len(cut_positions[node])):
				last_cut = 0
				if i > 0: last_cut = cut_positions[node][i-1]
				f.write(node + "_cut" + str(i) + "\t" + ">" + node + ":" + str(last_cut) + ":" + str(nodelens[node] - cut_positions[node][i]) + "\n")
			f.write(node + "_cut" + str(len(cut_positions[node])) + "\t" + ">" + node + ":" + str(cut_positions[node][-1]) + ":" + "0" + "\n")

with open(in_graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			if parts[1] not in cut_positions or len(cut_positions[parts[1]]) == 0:
				print(l.strip())
			else:
				for i in range(0, len(cut_positions[parts[1]])):
					last_cut = 0
					if i > 0: last_cut = cut_positions[parts[1]][i-1]
					node_name = parts[1] + "_cut" + str(i)
					node_seq = parts[2][last_cut:-(nodelens[parts[1]] - cut_positions[parts[1]][i])]
					print("S\t" + node_name + "\t" + node_seq + "\t" + "\t".join(parts[3:]))
				node_name = parts[1] + "_cut" + str(len(cut_positions[parts[1]]))
				node_seq = parts[2][cut_positions[parts[1]][-1]:]
				print("S\t" + node_name + "\t" + node_seq + "\t" + "\t".join(parts[3:]))
		if parts[0] == "L":
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode[1:] in cut_positions and len(cut_positions[fromnode[1:]]) > 0:
				if fromnode[0] == "<":
					fromnode = fromnode + "_cut0"
				else:
					fromnode = fromnode + "_cut" + str(len(cut_positions[fromnode[1:]]))
			if tonode[1:] in cut_positions and len(cut_positions[tonode[1:]]) > 0:
				if tonode[0] == ">":
					tonode = tonode + "_cut0"
				else:
					tonode = tonode + "_cut" + str(len(cut_positions[tonode[1:]]))
			parts[1] = fromnode[1:]
			parts[3] = tonode[1:]
			print("\t".join(parts))

