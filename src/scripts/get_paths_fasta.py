#!/usr/bin/env python3

import sys
import os
import graph_functions 
from scaffolding import path_storage


def main():
    rukki_path = sys.argv[1]
    fasta_file = sys.argv[2]
    graph = sys.argv[3]
    out_fasta = sys.argv[4]
    G = graph_functions.load_direct_graph(graph)
    path_storage = path_storage.PathStorage(G)
    path_storage.readFromFile(rukki_path)
    path_storage.writeFasta(fasta_file, out_fasta)