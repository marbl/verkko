#!/usr/bin/env python

import sys

graph_file = sys.argv[1]
node_lens_file = sys.argv[2]
rdna_nodes_file = sys.argv[3]
telomere_locations_file = sys.argv[4]
node_to_contig_file = sys.argv[5]

# graph to stdout

def revnode(n):
    assert len(n) >= 2
    assert n[0] == ">" or n[0] == "<"
    return (">" if n[0] == "<" else "<") + n[1:]

def canontip(left, right):
    fwstr = left + right
    bwstr = right + left
    if bwstr < fwstr: return (right, left)
    return (left, right)

def revcomp(s):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return "".join(comp[c] for c in s[::-1])

rdna = set()
with open(rdna_nodes_file) as f:
    for l in f:
        parts = l.strip().split('\t')
        rdna.add(parts[0])

translate_contig = {}
with open(node_to_contig_file) as f:
    for l in f:
        parts=l.strip().split(' ')
        if parts[0] == "path":
            translate_contig[parts[1]] = parts[2]

node_seqs = set()
telnodes = set()
node_lens = {}
edge_overlaps = {}

with open(node_lens_file) as f:
    for l in f:
        parts = l.strip().split('\t')
        node_lens[parts[0]] = int(parts[1])
with open(graph_file) as f:
    for l in f:
        parts = l.strip().split('\t')
        if parts[0] == "S":
            node_seqs.add(">" + parts[1])
            node_seqs.add("<" + parts[1])
        if parts[0] == 'L':
            fromnode = (">" if parts[2] == "-" else "<") + parts[1]
            tonode = ("<" if parts[4] == "-" else ">") + parts[3]
            if fromnode not in edge_overlaps:
                edge_overlaps[fromnode] = set() 
            edge_overlaps[fromnode].add(tonode)

with open(telomere_locations_file) as f:
    for l in f:
        parts = l.strip().split('\t')
        fromnode = ""

        #assert(parts[0] in translate_contig)
        if parts[0] not in translate_contig: continue
        graph_node = translate_contig[parts[0]]
        if int(parts[1]) < 10000:
            fromnode=">" + graph_node
            tonode=">telomere_"+graph_node+"_start"
        elif int(parts[2])+10000 > node_lens[parts[0]]:
            fromnode="<" + graph_node
            tonode=">telomere_"+graph_node+"_end"

        if fromnode == "": continue
        if tonode in telnodes: continue
        if fromnode in edge_overlaps:
            sys.stderr.write("Warning, didn't expect to have an overlap on telomere end of node %s line %s and fromnode is %s\n"%(graph_node, l.strip(), fromnode))
        telnodes.add(tonode)
        edge_overlaps[fromnode] = set()
        edge_overlaps[fromnode].add(tonode)
        edge_overlaps[tonode] = set()
        edge_overlaps[tonode].add(fromnode)

print("S\trDNA\t*\tLN:i:45000")
# output new telomere nodes
for t in telnodes:
    print("S\t%s\t*\tLN:i:6"%(t[1:]))
    if t in edge_overlaps:
        for l in edge_overlaps[t]:
            if l[0] == '>': ori="+"
            elif l[0] == '<':ori = "-" 
            print("L\t%s\t%s\t%s\t%s\t0M"%(t[1:], '+', l[1:], ori))

with open(graph_file) as f:
    for l in f:
        parts = l.strip().split('\t')
        if parts[0] == "S":
            if parts[1] not in rdna:
                print(l.strip())
        if parts[0] == 'L':
            if parts[1] not in rdna and parts[3] not in rdna:
                print(l.strip())
            elif parts[1] not in rdna and parts[3] in rdna:
                print("L\t%s\t%s\trDNA\t+\t0M"%(parts[1], parts[2]))
            elif parts[1] in rdna and parts[3] not in rdna:
                print("L\trDNA\t+\t%s\t%s\t0M"%(parts[3], parts[4]))
