#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions


#TODO: All classes and functions moved to verkko rep and should not be edited here but imported instead


#path nodes stored with orientation, "utig4-234+" tsv format
class pathStorage:
    #from IDs to array of oriented+- nodes
    paths = {}
    #ignoring overlaps but who cares
    path_lengths = {}
    def getPathById(self, path_id):
        return self.paths[path_id]
    
    def getPathIds(self):
        return self.paths.keys()

    def getEdgeSequenceById(self, path_id):
        res = []
        for edge in self.paths[path_id]:
            #removing +- and gaps
            if edge[0] != "N":
                res.append(edge[:-1])
        return res

    def addPath(self, line, G):
#        print (line)
        arr = line.strip().split()
        if len(arr) < 3:
            print (f"wrong path {line}")
            exit()
        separators = ">|<|\[|\]|,"
        edges = re.split(separators, arr[1])
        total_l = 0
        for edge in edges:
            node = edge[:-1]
            if node in G.nodes:
                total_l += G.nodes[node]['length']
        edges = list(filter(None, edges))
#        print (edges)
        self.paths[arr[0]] = edges
        self.path_lengths[arr[0]] = total_l



def get_rdna_mashmap(hicrun_dir):
    IDY_THRESHOLD = 0.95
    
    hic_stage = hicrun_dir

    if not os.path.exists(os.path.join(hic_stage, "rdna_mashmap.out")):
        sys.stderr.write(os.path.join(hic_stage, "rdna_mashmap.out") + " WTF") 
        sys.stderr.flush()
            #TODO This block should never be run in verkko pipeline
        script_file = os.path.join(hic_stage, "rdna_mashmap.sh")    
        with open (script_file, 'w') as scr_file:
            scr_file.write("#!/usr/bin/bash\n")    
            scr_file.write ("module load mashmap\n")
            scr_file.write(f"mashmap -s 5000  --pi 95.0 -t 24 -r /data/antipovd2/devel/ribotin/template_seqs/chm13_rDNAs.fa -q {hic_stage}/unitigs.fasta -o {hic_stage}/rdna_mashmap.out")
            scr_file.close()
        os.system(f"bash {script_file}")
    
    rdna_nodes = set()
    node_cov = {}
    node_len = {}
    for l in open(os.path.join(hic_stage, "rdna_mashmap.out")):
        arr = l.strip().split()        
        node = arr[0]
        if len(arr) > 10:
            idy = float(arr[12].split(":")[2])
        else:
            idy = float(arr[9])/100
        if not node in node_cov:
            node_cov[node] = []
        if idy < IDY_THRESHOLD:
            continue 
        node_len[node] = int(arr[1])
        node_cov[node].append([int (arr[2]), 1])
        node_cov[node].append([int (arr[3]), -1])                        
    rdna_nodes = set()
    for node in node_len.keys():
        covered = 0
        prev = 0
        state = 0
        pos = sorted(node_cov[node], key=lambda x: x[0])
        
        for p in pos:
            if state > 0:
                covered += p[0] - prev
            prev = p[0]
            state += p[1]
        if covered > 0.1 * node_len[node]:
            rdna_nodes.add(node)
        else:
            print (f"node {node} skipping. coverage {covered} of {node_len[node]}")
    return rdna_nodes

def rc_path(path):
    res = []
    path.reverse()
    for node in path:
        if node[-1] == "+":
            res.append(node[:-1]+"-")
        elif node[-1] == "-":
            res.append(node[:-1]+"+")   
        else:
            res.append(node)
    return res
    
#Should not be run in normal verkko pipeline
def get_rdna_ribotin(hicrun_dir):
    rdna_file = os.path.join(hicrun_dir, os.pardir, "rdnanodes.txt")     
    if not(os.path.exists(rdna_file)):
        hpc_graph = os.path.join(hicrun_dir, os.pardir, 'assembly.homopolymer-compressed.gfa')
        os.system(f"cp  {os.path.join(hicrun_dir, os.pardir, '5-untip/unitig-unrolled-unitig-unrolled-popped-unitig-normal-connected-tip.gfa')} {hpc_graph}")
        os.system(f"/data/antipovd2/devel/ribotin/bin/ribotin-verkko  -i {os.path.join(hicrun_dir, os.pardir)} -o {os.path.join(hicrun_dir, os.pardir,'ribo_test')} -x human --do-ul no -t 20")
        os.system(f"cat {os.path.join(hicrun_dir, os.pardir)+'/ribo_test*/nodes.txt'} > {rdna_file}")
    rdna_nodes = set()    
    for l in open(rdna_file):
        arr = l.strip().split()
        for e in arr:
            rdna_nodes.add(e)
    return rdna_nodes

def get_telomeric_nodes(telomere_locations_file, G):
    aux_tel_nodes = set()
    real_tel_nodes = set()
    with open(telomere_locations_file) as f:
        for l in f:
            parts = l.strip().split('\t')
            telnode = ""
            graph_node = parts[0]
            if int(parts[1]) < 10000:
                telnode=">telomere_"+graph_node+"_start"
                G.add_node(telnode, length=0, coverage=0)
                G.add_edge(telnode, graph_node)
            elif int(parts[2])+10000 > G.nodes[graph_node]['length']:
                telnode=">telomere_"+graph_node+"_end"
                G.add_node(telnode, length=0, coverage=0)
                G.add_edge(graph_node, telnode)            
            aux_tel_nodes.add(telnode)
            real_tel_nodes.add(graph_node)
    return aux_tel_nodes, real_tel_nodes

def read_rukki_paths(rukki_file, G):
    res = pathStorage()
    for line in open (rukki_file):
        arr = line.strip().split()
        if arr[0] == "name":
            continue
        res.addPath(line.strip(), G)
    return res

#
#Map from path-id to direction (suspected telomere left, rdna right - true - otherwise false) 
def get_arm_paths_ids(telomere_locations_file, G, paths, rdna_nodes, dists, min_path_len, max_path_len, ignore_telomeres):
    arm_paths_ids = {}
    aux_tel_nodes, real_tel_nodes = get_telomeric_nodes(telomere_locations_file, G)
    for id in paths.paths:
        total_l = 0
        for node in paths.getEdgeSequenceById(id):
#        print (paths.getPathById(id))
            if node in G.nodes:
                total_l += G.nodes[node]['length']
#short arm length condition
    #should we take in account already resolved rDNAs??
        if total_l < min_path_len or total_l > max_path_len:
            continue
        path = paths.getEdgeSequenceById(id)
        close_telo = False
        pathends = {path[0], path[-1]}
        for pathend in pathends:
#TODO better use of coords            
            if pathend in real_tel_nodes:
                print (id + " in telomeres check")
                close_telo = True
#oriented ones, so telnodes
            for tels in aux_tel_nodes:
                if tels in dists and pathend in dists[tels] and dists[tels][pathend] < G.nodes[pathend]['length'] + 200000:
                    print (id + " near telomeres check")
                    close_telo = True
                    continue
        close_rdna = False
#TODO: should this be here for long arms?        
        for node in path:
            if node in rdna_nodes:
                print (id + " rDNA in path")
                close_rdna = True

        direction = True
        closest_rdna_node = ""
        closest_rdna_dist = 999999999
        closest_suff_node = ""
        for suff in pathends:
            for t in rdna_nodes:
                if t in dists and suff in dists[t]:
                    if dists[t][suff] < closest_rdna_dist:
                        closest_rdna_dist = dists[t][suff]
                        closest_rdna_node = t                        
                        closest_suff_node = suff
                    if dists[t][suff]< G.nodes[suff]['length']+ G.nodes[t]['length'] + 2000000:
                        close_rdna = True
                        if suff == path[0]:
                            direction = False
                        print (id + " rDNA near path")
                        break
        if not close_rdna and closest_rdna_dist < 999999999:
            print (id + " rDNA not so near path, closest node " + closest_rdna_node + " dist " + str(closest_rdna_dist) + " total_l " +str(G.nodes[closest_suff_node]['length']+ G.nodes[closest_rdna_node]['length']))
        if close_rdna and (close_telo or ignore_telomeres):
            arm_paths_ids[id] = direction
    return arm_paths_ids


def score(node, path, multiplicities, weights_map, use_multiple):
    res = 0
    for edge_ext in path:
        edge = edge_ext
        if multiplicities[edge] == 1 or use_multiple:
            if edge in weights_map.keys():
                if node in weights_map[edge].keys():                
                    res += weights_map[edge][node]
#        else:
#            print (f"multiplicity {multiplicities[edge]} for {edge} skipped")
    return res

#estimating median background links count for short arm, assuming that majority of long nodes have no connection to that short arm
def get_background_score(short_path, weights_map, G):    
    MIN_LONG_LEN = 10000000
    avg_links = []
    for node in G.nodes:
        if G.nodes[node]['length'] > MIN_LONG_LEN and not node in short_path and node in weights_map.keys():
            sum_link = 0
            sum_len = 0
            for short_node in short_path:
                if short_node in weights_map[node]:
                    sum_link += weights_map[node][short_node]
            avg_links.append(sum_link / G.nodes[node]['length'])
    print (f"path {short_path} avg links {len(avg_links)} {avg_links[int(len(avg_links)/2)]}")
    if len(avg_links) == 0:
        return 0
    avg_links.sort()
    return avg_links[int(len(avg_links)/2)]

#for components we allow to use nodes present in multiple paths, for haplo paths we do not
def scorepath(short_arm_path, long_path, mults, weights_map, use_multiple):    
    res = 0
    for node in short_arm_path:
        if mults[node] == 1:
            res += score(node, long_path, mults, weights_map, use_multiple)
    return res

def path_length(path, G):
    res = 0
    for edge in path:
        res += G.nodes[edge]['length']
    return res

def get_best_path(short_id, path_container, paths, longarm_to_component, multiplicities, weights_map, G):
    path_scores = {}
    best_arms = []
    for long_arm in path_container.keys():
        path_scores[long_arm] = scorepath(paths.getEdgeSequenceById(short_id), path_container[long_arm], multiplicities, weights_map, True)
        best_arms.append(long_arm)
    best_arms.sort(key=lambda x: path_scores[x], reverse=True)
    if len(best_arms) == 1:
        return best_arms[0]
    background_score = get_background_score(paths.getEdgeSequenceById(short_id), weights_map, G)
    if path_scores[best_arms[0]] - (path_length(path_container[best_arms[0]], G) * background_score) < 2 * (path_scores[best_arms[1]] - path_length(path_container[best_arms[1]], G) * background_score):
        print (f"No evident majority, best options are {best_arms[0]} and {best_arms[1]}, checking further")
        if longarm_to_component[best_arms[0]] != longarm_to_component[best_arms[1]]:
            #not homologous chromosomes, we can't resolve it
            return "Unclear"
        het_path_scores = {}
        for long_arm in path_container.keys():
            het_path_scores[long_arm] = scorepath(paths.getEdgeSequenceById(short_id), path_container[long_arm], multiplicities, weights_map, False)
        if het_path_scores[best_arms[0]] - (path_length(path_container[best_arms[0]], G) * background_score) < 2 * (het_path_scores[best_arms[1]] - path_length(path_container[best_arms[1]], G) * background_score):
            print(f"still no evident majority, best options are {best_arms[0]} and {best_arms[1]}, scores {path_scores[best_arms[0]]} and {path_scores[best_arms[1]]}, het scores {het_path_scores[best_arms[0]]} and {het_path_scores[best_arms[1]]}")
            return "Unclear"
    return best_arms[0]
