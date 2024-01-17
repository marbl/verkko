#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions as gf 
import rdna_scaff_functions as sf

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
if os.path.exists(old_rukki_tsv_file):
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

#read graph
#add_telomeric_nodes
#construct new graph (vertex = start/end of edge), ignoring overlaps
#dijkstra 
#condition 

#not expexted to be run in main pipeline
if debug and os.path.exists(rdna_file): 
    print("Somehow you are deeply in verkko's debug, is it intentional?")
    rdna_nodes_mashmap = sf.get_rdna_ribotin(hicrun_dir)
    rdna_nodes2= sf.get_rdna_mashmap(hicrun_dir)
    print ("Ribo non mashmap")
    print(rdna_nodes_mashmap - rdna_nodes2)
    print ("Mashmap non ribo")
    print(rdna_nodes2 - rdna_nodes_mashmap)

rdna_nodes = sf.get_rdna_mashmap(hicrun_dir)


G = nx.DiGraph()
#graph_f fails?
#TODO should use graph_functions!
gf.load_direct_graph(gfa_file, G)
#double distance between middle of nodes
dists = dict(nx.all_pairs_dijkstra_path_length(G, weight='mid_length'))
paths = sf.read_rukki_paths(rukki_tsv_file, G)

short_arm_paths_ids = sf.get_arm_paths_ids(telomere_locations_file, G, paths, rdna_nodes, 100000, 5000000, False)
experimental_long_arm_path_ids = sf.get_arm_paths_ids(telomere_locations_file, G, paths, rdna_nodes, 10000000, 999999999, True)
same_component_long_arm_path_ids = sf.get_same_component_paths(short_arm_paths_ids, G, paths, 10000000,999999999)
for idd in same_component_long_arm_path_ids.keys():
    if idd not in experimental_long_arm_path_ids:
        experimental_long_arm_path_ids[idd] = same_component_long_arm_path_ids[idd]
        print (f"Component-based added path {idd}")
print (f"short arms PATH number{len(short_arm_paths_ids)}")
for idd in short_arm_paths_ids.keys():
    print(paths.getPathById(idd))

print (short_arm_paths_ids)
#parse translation_hap files

#getting long arms
#haplotype2-0000054      chr22   45121307        51324926
#chr names in general case not available

#Should not be run in normal runs
if debug and os.path.exists(translation_hap1_file) and os.path.exists(translation_paths_file):
    print("Somehow you are deeply in verkko's debug, is it intentional?")
    long_arm_paths = {}
    multiplicities = {}
    acros = {"chr13", "chr14", "chr15", "chr21", "chr22"}
    long_arms = {}
    for i in range (1, 3):
        trans_file = os.path.join (hicrun_dir, os.pardir, "translation_hap" + str(i))
        for line in open(trans_file):
            [id, chr_num, plen, chr_len] = line.split()
            if chr_num in acros and int(plen) < int(chr_len):
                long_arm = chr_num + "_" + str(i)
                #let's use longest assignment that is shorter than fixed known chr length
                if long_arm in long_arms and long_arms[long_arm][1]> int(plen):
                    continue
                long_arms[long_arm] =[id, int(plen)]
    chr_to_rukkiid = {}
    rukkiid_to_name = {}
    for line in open(translation_paths_file):
        arr = line.strip().split()
        if len(arr) == 3 and arr[0] == "path":
            chr_to_rukkiid[arr[1]] = arr[2]     
            rukkiid_to_name[arr[2]] = arr[1]

    print ("long arm rukki ids:")
    mapping_based_long_arm_path_ids = set()
    for long_arm in sorted(long_arms.keys()):
        mapping_based_long_arm_path_ids.add(chr_to_rukkiid[long_arms[long_arm][0]])
    print ("experimental but not mapping based:")
    print (experimental_long_arm_path_ids.keys() - mapping_based_long_arm_path_ids)

    print ("Mapping based but not experimental based:")
    print (mapping_based_long_arm_path_ids - experimental_long_arm_path_ids.keys())


long_arm_paths = {}
for long_arm in sorted(experimental_long_arm_path_ids.keys()):
    long_arm_paths[long_arm] = paths.getPathById(long_arm)
print (f"Total long arms found {len(experimental_long_arm_path_ids.keys())}")


multiplicities = {}
for path_id in paths.getPathIds():
    for edge in paths.getEdgeSequenceById(path_id):
        for dir in ["+", "-"]:
            node = edge + dir      
            if not node in multiplicities:
                multiplicities[node] = 0
            multiplicities[node] += 1



components = {}
component_id = 0
node_to_component = {}
for line in open(hicverkko_log):
    arr = line.strip().split()
    if arr[0] == "Connected":
        component_id += 1
        utigs = re.findall(r"'(utig\d+-\d+)'", line)
        for utig in utigs:
            node_to_component[utig] = component_id
        components[component_id] = utigs

longarm_to_component = {}
longarm_component_nodes = {}
component_to_longarm = {}
for long_arm in sorted(experimental_long_arm_path_ids.keys()):
    path = long_arm_paths[long_arm]
    print (f"long arm {long_arm} path {path}")
    comp_w = {}
    for component_id in components.keys():
        comp_w[component_id] = 0
    for node in path:
        if node in node_to_component:
            comp_w[node_to_component[node]] += G.nodes[node+'+']['length']
    majority_component = max(comp_w, key=comp_w.get)
    longarm_to_component[long_arm] = majority_component
    if not majority_component in component_to_longarm:
        component_to_longarm[majority_component] = []
        longarm_component_nodes[majority_component] = []
    component_to_longarm[majority_component].append(long_arm)
    longarm_component_nodes[majority_component].extend(long_arm_paths[long_arm])

print("long arm components")
print(longarm_to_component)    


weights_map = {}
for line in open (hic_byread):
    arr = line.strip().split()
    if not arr[1]+'+' in weights_map:
        for dir in ["+", "-"]:
            weights_map[arr[1] + dir] = {}
    if not arr[2]+'+' in weights_map:
        for dir in ["+", "-"]:
            weights_map[arr[2] + dir] = {}
    for dir1 in ["+", "-"]:
        for dir2 in ["+", "-"]:
            weights_map[arr[1]+dir1][arr[2] + dir2] =  int(arr[3])
            weights_map[arr[2]+dir1][arr[1] +dir2] =  int(arr[3])
         
#path scoring
final_assignment = {}

tsv_out = {}
#if one long assigned to multiple short - do not change any and report weirdness
long_arm_multiplicities = {}
for short_id in short_arm_paths_ids.keys():
    short_arm_nodes = paths.getPathById(short_id)
    #    print (weights_map[short_id])
    max_score = 0
    max_long_arm = "none"
    best_path = sf.get_best_path(short_id, long_arm_paths, paths, longarm_to_component, multiplicities, weights_map, G)
    final_assignment[short_id] = best_path
    if not best_path in long_arm_multiplicities:
        long_arm_multiplicities[best_path] = 0
    long_arm_multiplicities[best_path] += 1
for short_id in final_assignment.keys():
    best_path = final_assignment[short_id]
    if best_path != "Unclear":
        if long_arm_multiplicities[best_path] > 1:
            error_str = f"Multiple short arms assigned to one long arm {best_path}, not scaffolding!"
            for id in final_assignment.keys():
                if final_assignment[id] == best_path:
                    error_str += f" {id}"
                    final_assignment[id] = "Unclear"
            print (error_str)
        else:
            to_output = []
            if experimental_long_arm_path_ids[best_path]:
                to_output.extend(paths.getPathById(best_path))
            else:
                to_output.extend(sf.rc_path(paths.getPathById(best_path)))
            to_output.append("[N1000001N:rdna]")
            if short_arm_paths_ids[short_id]:
                to_output.extend(sf.rc_path(paths.getPathById(short_id)))
            else:
                to_output.extend(paths.getPathById(short_id))
            res_path = ','.join(to_output)
            tsv_out[best_path] = res_path
            tsv_out[short_id] = "NONE"

with open(scaff_rukki_tsv_file, "w") as out_tsv:        
    for line in open(rukki_tsv_file):
        arr = line.strip().split()
        if arr[0] in tsv_out:
            if tsv_out[arr[0]] == "NONE":
                continue
            else: 
                out_tsv.write(arr[0] + "\t" + tsv_out[arr[0]] + "\t" + arr[2] + "\n")
        else:
            out_tsv.write(line)
with open(scaff_rukki_gaf_file, "w") as out_gaf:        
    for line in open(scaff_rukki_tsv_file):
        arr = line.strip().split()
        out_gaf.write(arr[0] + "\t" + gf.tsv2gaf(arr[1]) + "\t" + arr[2] + "\n")
        
#Additional output for debugging only
scores = {}
assgnd = 0
for short_id in short_arm_paths_ids.keys():
    print (f"Assigned PATH: {short_id} --- {final_assignment[short_id]}")
    if final_assignment[short_id] != "Unclear":
       assgnd += 1
    scores[short_id] = {}
    for long_arm in long_arm_paths.keys():
        scores[short_id][long_arm] = sf.scorepath(paths.getPathById(short_id), long_arm_paths[long_arm], multiplicities, weights_map, True)
    print (f"All scores path using homozygous: {short_id} --- {scores[short_id]}\n")
    
    scores[short_id] = {}
    for long_arm in long_arm_paths.keys():
        scores[short_id][long_arm] = sf.scorepath(paths.getPathById(short_id), long_arm_paths[long_arm], multiplicities, weights_map, False)
    print (f"All scores NOT using homozygous: {short_id} --- {scores[short_id]}\n")
        
print (f"total paths assigned {assgnd} of {len(short_arm_paths_ids)}")
