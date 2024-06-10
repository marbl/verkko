#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions as gf 
import logging
from scaffolding import scaffold_graph, logger_wrap, path_storage

#get all rdna nodes from mashmap
#get nodes from seqtk telo 
#select nodes that are between 500K and 2M and close to both rDNA and telomeres to be short arm

debug = False

run_dir = sys.argv[1]

#TODO: deprectated
telomere_script = "/data/antipovd2/devel/scripts/verkko_output_analysis//getOnlyTelomere.sh"
hap_assign_script = "/data/antipovd2/devel/scripts/verkko_output_analysis//get.sh"

hicrun_dir = os.path.realpath(run_dir)

telomere_locations_file = os.path.join(hicrun_dir,  "unitigs.telo")
uncompressed_nodes = os.path.join(hicrun_dir,  "unitigs.fasta")
#node_to_contig_file = os.path.join(oldrun_nodes_dir, "6-layoutContigs", "unitig-popped.layout.scfmap")

gfa_file = os.path.join(hicrun_dir, "unitigs.hpc.noseq.gfa")
rukki_tsv_file = os.path.join(hicrun_dir,  "rukki.paths.tsv")
old_rukki_tsv_file = os.path.join(hicrun_dir,  "prescaf_rukki.paths.tsv")
if debug and os.path.exists(old_rukki_tsv_file):
    os.system(f"cp {old_rukki_tsv_file} {rukki_tsv_file}")
    os.system (f'cp {os.path.join(hicrun_dir,  "prescaf_rukki.paths.gaf")} {os.path.join(hicrun_dir,  "rukki.paths.gaf")}')
    print ("Scaffolder rerunning, restoring prescaf rukki file")
scaff_rukki_tsv_file = os.path.join(hicrun_dir,  "scaff_rukki.paths.tsv")
scaff_rukki_gaf_file = os.path.join(hicrun_dir,  "scaff_rukki.paths.gaf")
translation_paths_file = os.path.join (hicrun_dir, "final_contigs/6-layoutContigs/unitig-popped.layout.scfmap")
translation_hap1_file = os.path.join (hicrun_dir, os.pardir, "translation_hap1")
rdna_file = os.path.join(hicrun_dir, os.pardir, "rdnanodes.txt")     

hic_byread = os.path.join(hicrun_dir, "hic.byread.compressed")

hicverkko_log = os.path.join(hicrun_dir, "hicverkko.log")


G = nx.DiGraph()
gf.load_direct_graph(gfa_file, G)
logger = logger_wrap.initLogger(os.path.join(hicrun_dir, "scaffolding.log"))
path_mashmap = os.path.join(hicrun_dir, "paths2ref.mashmap")
rukki_path_storage = path_storage.PathStorage(G)
rukki_path_storage.readFromFile(rukki_tsv_file)
#sf.try_to_scaff(rukki_paths, telomere_locations_file, os.path.join(hicrun_dir, "hic_mapping.byread.output"), os.path.join(hicrun_dir, "unitigs.matches"), G, indirectG, uncompressed_nodes)
sg = scaffold_graph.ScaffoldGraph(rukki_path_storage, telomere_locations_file, os.path.join(hicrun_dir, "hic_mapping.byread.output"), os.path.join(hicrun_dir, "unitigs_nonhpc50.mashmap"), G, uncompressed_nodes, path_mashmap, logger) 
res = sg.generateScaffolds()

