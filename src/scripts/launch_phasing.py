#!/usr/bin/env python3
import shutil
import sys
import os
import random
from os import listdir
from os.path import isfile, join
import cluster
if __name__ == "__main__":
    if len(sys.argv) < 4:
        print(f'Usage: {sys.argv[0]} <no_rdna_tangle> <uneven_depth> <output_dir>')
        print(f'no_rdna & uneven_depth are bool variables, determing whether to switch off rdna_tangle based heuristics and whether coverage is highly uneven. Defaults: False, False')
        exit()
    #run clustering
    #python3 cluster.py unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.noseq.gfa unitig-popped-unitig-normal-connected-tip.homopolymer-compressed.matches hic_mapping.byread.output > cluster.out 2> cluster.err

    #convert to rukki inputs
    #sh convert.sh

    # run rukki
    #sh rukki.sh
    cur_dir = os.path.abspath(os.path.dirname(__file__))
    no_rdna = eval(sys.argv[1])
    uneven_depth = eval(sys.argv[2])
    output_dir = sys.argv[3]

    os.makedirs(output_dir, exist_ok=True)
    matches_file = os.path.join(output_dir, "unitigs_hpc50.mashmap")
    hic_file = os.path.join(output_dir, "hic_mapping.byread.output")
    if not os.path.exists(hic_file):
        hic_file = os.path.join(output_dir, "hic.byread.compressed")
    compressed_hic = os.path.join(output_dir, "hic.byread.compressed")
    if os.path.exists(compressed_hic):
        hic_file = compressed_hic

    noseq_gfa = os.path.join(output_dir, "unitigs.hpc.noseq.gfa")
    cluster.run_clustering(noseq_gfa, matches_file, hic_file, output_dir, no_rdna, uneven_depth)
#Saved for further debug


#    tst_dir = os.path.join(output_dir, "tst_run")   
#    os.makedirs(tst_dir, exist_ok=True)
#    cluster.run_clustering(noseq_gfa, matches_file, hic_file, tst_dir, no_rdna, uneven_depth)
