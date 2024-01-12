#!/usr/bin/env python

import sys
import argparse

# Create the parser
parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument('-r', '--rdna', required=False, help='rDNA node list, from mash screen')
parser.add_argument('-t', '--telo', required=False, help='Telomere location in bed format, from seqtk telo -d 10000')
parser.add_argument('-g', '--gfa', default='assembly.homopolymer-compressed.noseq.gfa', help='assembly.homopolymer-compressed.noseq.gfa')
parser.add_argument('-s', '--scfmap', default='assembly.scfmap', help='assembly.scfmap')
parser.add_argument('-p', '--paths', default='assembly.paths.tsv', help='assembly.paths.tsv')
parser.add_argument('-o', '--output', default='assembly.homopolymer-compressed.noseq.telo_rdna.gfa', help='Modified GFA with added telomere nodes and an rDNA node for a mok-collapsed rDNA (default: assembly.homopolymer-compressed.noseq.telo_rdna.gfa)')
parser.add_argument('-c', '--colors', default='assembly.colors.telo_rdna.csv', help='Colors for the telomere and rDNA nodes to go along with the modified GFA. Nodes will be colored green for telomeres, purple for rDNA node. If assembly.colors.csv exist, it will be copied to this output with the new colors appended at the end. (default: assembly.colors.telo_rdna.csv)')

# Parse the arguments
args = parser.parse_args()

rdna_nodes_file = args.rdna     # rdna.screennodes.list from mash
telomere_file = args.telo       # telomere from seqtk telo -d 10000

graph_file = args.gfa           # assembly.homopolymer-compressed.noseq.gfa
node_to_path_file = args.scfmap # assembly.scfmap
path_to_nodes_file = args.paths # assembly.paths.tsv

# Check if either rdna or telo is provided
if not args.rdna and not args.telo:
    parser.error("--rdna or --telo must be provided.")

telo_rdna_gfa_file = open(args.output, "w")

telo_rdna_color_file = open(args.colors, "w")
padding="\t0\t0\t0:0\t"

try:
    with open("assembly.colors.csv", "r") as f:
        for l in f:
            telo_rdna_color_file.write(l)
except:
    sys.stderr.write("No assembly.colors.csv file found. Creating a simple colored csv file\n")
    padding = "\t"
    telo_rdna_color_file.write("node\tcolor\n")

telo_col = "#008000" # green
rdna_col = "#A020F0" # purple

rdna = set()
try:
    with open(rdna_nodes_file, "r") as f:
        for l in f:
            parts = l.strip().split('\t')
            rdna.add(parts[0])
except:
    sys.stderr.write("No rDNA nodes found. Skipping rDNA node creation\n")

has_rDNA = True if len(rdna) > 0 else False

# contig : path
contig_to_path = {}
translate_contig = {}
with open(node_to_path_file, "r") as f:
    for l in f:
        parts=l.strip().split(' ')
        if parts[0] == "path": 
            # parts[1] = contig name
            # parts[2] = path name
            contig_to_path[parts[1]] = parts[2]
            #for hi-c runs
            contig_to_path[parts[2]] = parts[2]



# path : nodes
path_to_nodes = {}
with open(path_to_nodes_file, "r") as f:
    for l in f:
        parts=l.strip().split('\t')
        # parts[0] = path name
        # parts[1] = nodes in path, "," separated
        path_to_nodes[parts[0]] = [ parts[1].split(',')[0], (parts[1].split(','))[-1] ]

# node_seqs = set()
telnodes = set()
edge_overlaps = {}

with open(graph_file, "r") as f:
    for l in f:
        parts = l.strip().split('\t')
        if parts[0] == 'L':
            fromnode = (">" if parts[2] == "-" else "<") + parts[1]
            tonode = ("<" if parts[4] == "-" else ">") + parts[3]
            if fromnode not in edge_overlaps:
                edge_overlaps[fromnode] = set() 
            edge_overlaps[fromnode].add(tonode)

with open(telomere_file, "r") as f:
    for l in f:
        parts = l.strip().split('\t')
        fromnode = ""

        #assert(parts[0] in translate_contig)
        if parts[0] not in contig_to_path:
            sys.stderr.write("Warning, no path available for contig %s line %s\n"%(parts[0], l.strip()))
            continue

        path = contig_to_path[parts[0]]
        #for hic
        if path.startswith("utig"):
            graph_node = [path, path+"+"]
        else:
            graph_node = path_to_nodes[path]
        if int(parts[1]) < 10000:
            telnode=graph_node[0]
            # reverse the beginning node direction
            fromnode= (">" if telnode[-1] == "+" else "<") + telnode[0:-1]
            tonode=">telomere_" + telnode + "_start"
            telo_rdna_color_file.write("%s%s%s\n"%(tonode[1:], padding, telo_col))
        elif int(parts[2])+10000 > int(parts[3]):
            telnode=graph_node[1]
            fromnode=(">" if telnode[-1] == "-" else "<") + telnode[0:-1]
            tonode=">telomere_" + telnode + "_end"
            telo_rdna_color_file.write("%s%s%s\n"%(tonode[1:], padding, telo_col))

        if fromnode == "": continue
        if tonode in telnodes: continue
        if fromnode in edge_overlaps:
            sys.stderr.write("Warning, didn't expect to have an overlap on telomere end of node %s line %s and fromnode is %s\n"%(graph_node, l.strip(), fromnode))
        telnodes.add(tonode)
        edge_overlaps[fromnode] = set()
        edge_overlaps[fromnode].add(tonode)
        edge_overlaps[tonode] = set()
        edge_overlaps[tonode].add(fromnode)

if (has_rDNA):
    telo_rdna_gfa_file.write("S\trDNA\t*\tLN:i:45000\n")
    telo_rdna_color_file.write("rDNA%s%s\n"%(padding, rdna_col))

# output new telomere nodes
for t in telnodes:
    telo_rdna_gfa_file.write("S\t%s\t*\tLN:i:6\n"%(t[1:]))
    if t in edge_overlaps:
        for l in edge_overlaps[t]:
            if l[0] == '>': ori="+"
            elif l[0] == '<':ori = "-" 
            telo_rdna_gfa_file.write("L\t%s\t%s\t%s\t%s\t0M\n"%(t[1:], '+', l[1:], ori))

with open(graph_file, "r") as f:
    for l in f:
        parts = l.strip().split('\t')
        if parts[0] == "S":
            if parts[1] not in rdna:
                telo_rdna_gfa_file.write(l.strip()+"\n")
        if parts[0] == 'L':
            if parts[1] not in rdna and parts[3] not in rdna:
                telo_rdna_gfa_file.write(l.strip()+"\n")
            elif parts[1] not in rdna and parts[3] in rdna:
                telo_rdna_gfa_file.write("L\t%s\t%s\trDNA\t+\t0M\n"%(parts[1], parts[2]))
            elif parts[1] in rdna and parts[3] not in rdna:
                telo_rdna_gfa_file.write("L\trDNA\t+\t%s\t%s\t0M\n"%(parts[3], parts[4]))
