#!/usr/bin/env python3

import sys
import os
import graph_functions as gf
from scaffolding import path_storage
import networkx as nx

rukki_path = sys.argv[1]
fasta_file = sys.argv[2]
graph_file = sys.argv[3]
out_fasta = sys.argv[4]
out_id_list = sys.argv[5]
path_len_cutoff = int(sys.argv[6])
path_end_cutoff = int(sys.argv[7])
fasta_uncompressed = sys.argv[8]
G = nx.DiGraph()
gf.load_direct_graph(graph_file, G)
path_storage = path_storage.PathStorage(G)
path_storage.readFromFile(rukki_path)
path_storage.writePathAsFasta(fasta_file, out_fasta)
interesting = set()
lens = gf.get_lengths(fasta_uncompressed)

#TODO short ignored node disappeared here, possibly should be returned sometimes
IGNORED_LEN = 0

multiplicities = path_storage.getEdgeMultiplicities()
for path_id in path_storage.getPathIds():            
    #TODO possibly reconsider this condition; /2 because of compressed/uncompressed
    if path_storage.getLength(path_id) < path_len_cutoff:
        continue
    total_len = 0
    for or_node in path_storage.getPathById(path_id):
        nor_node = or_node.strip('-+')
        
        if nor_node in lens and lens[nor_node] > IGNORED_LEN and or_node in multiplicities and multiplicities[or_node] == 1:
            total_len += lens[nor_node]
    before = 0
    after = total_len
    for or_node in path_storage.getPathById(path_id):        
        nor_node = or_node.strip('-+')
        if nor_node in lens and lens[nor_node] > IGNORED_LEN and or_node in multiplicities and multiplicities[or_node] == 1:
            after -= lens[nor_node]
            if before < path_end_cutoff or after < path_end_cutoff:
                interesting.add(nor_node)
            before += lens[nor_node]

with open(out_id_list, 'w') as id_list_file:
    for node in sorted(interesting):
        id_list_file.write(node + "\n")
