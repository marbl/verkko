#!/usr/bin/env python3

import sys
import os
import graph_functions 
from scaffolding import path_storage
import networkx as nx

rukki_path = sys.argv[1]
fasta_file = sys.argv[2]
graph_file = sys.argv[3]
out_fasta = sys.argv[4]
G = nx.DiGraph()
graph_functions.load_direct_graph(graph_file, G )
path_storage = path_storage.PathStorage(G)
path_storage.readFromFile(rukki_path)
path_storage.writePathAsFasta(fasta_file, out_fasta)
