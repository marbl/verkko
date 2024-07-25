#!/usr/bin/env python

import sys
import graph_functions as gf

graph_file = sys.argv[1]
# graph to stdout

not_tips = set()
node_seqs = set()
edge_overlaps = {}

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			node_seqs.add(">" + parts[1])
			node_seqs.add("<" + parts[1])
		if parts[0] == 'L':
			fromnode = (">" if parts[2] == "+" else "<") + parts[1]
			tonode = ("<" if parts[4] == "+" else ">") + parts[3]
			if fromnode not in edge_overlaps:
				edge_overlaps[fromnode] = set() 
			edge_overlaps[fromnode].add(tonode)
			if parts[1] == parts[3]: # we will skip self-edges for this consideration, if you have a loop at a gap, try to patch it anyway
				continue
			not_tips.add(fromnode)
			not_tips.add(tonode)

tips=node_seqs.difference(not_tips)
toadd=set()

with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == 'L':
			fromnode = ("<" if parts[2] == "+" else ">") + parts[1]
			tonode = (">" if parts[4] == "+" else "<") + parts[3]
			if fromnode in tips:
				tonode = gf.revnode(tonode)
				sys.stderr.write("The node %s is a tip and the list of nodes matching %s is %s\n"%(fromnode, tonode, edge_overlaps[tonode]))
				for i in edge_overlaps[tonode]:
					i=gf.revnode(i)
					if i not in edge_overlaps: continue
				toadd.update(edge_overlaps[tonode])
			elif tonode in tips:
				fromnode = gf.revnode(fromnode)
				sys.stderr.write("The node %s is a tip and the list of nodes matching %s is %s\n"%(tonode, fromnode, edge_overlaps[fromnode]))
				for i in edge_overlaps[fromnode]:
					i=gf.revnode(i)
					if i not in edge_overlaps: continue
				toadd.update(edge_overlaps[fromnode])

tips.update(toadd)
with open(graph_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] == "S":
			if ">" + parts[1] in tips or "<" + parts[1] in tips:
				print(">%s\n%s"%(parts[1], parts[2]))
