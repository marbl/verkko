#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions
import logging
from scaffolding import path_storage


    
#Actyally should be processed the same way as telos, with fake nodes of zero length
#This will reduce code duplication and will allow to use "one-sided" rdna nodes that may be useful
def get_rdna_mashmap(hicrun_dir):
    IDY_THRESHOLD = 0.95
    
    hic_stage = hicrun_dir

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
            for direction in ["+", "-"]:
                rdna_nodes.add(node+ direction)
        else:
            print (f"node {node} skipping. coverage {covered} of {node_len[node]}")
    return rdna_nodes

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
            for direction in ["+", "-"]:
                rdna_nodes.add(e + direction)
    return rdna_nodes

#return paths that are reachable from any of the rdna nodes and within length limits
def get_same_component_paths(short_arm_path_ids, G, paths, min_path_len, max_path_len):

    res = {}
    dists = dict(nx.all_pairs_dijkstra_path_length(G, weight='mid_length'))
    #nodes near short arm telomeres
    short_arm_starts = []
    for id in short_arm_path_ids:
        if short_arm_path_ids[id]:
            short_arm_starts.append(paths.getPathById(id)[0])
        else:   
            short_arm_starts.append(rc_path(paths.getPathById(id))[0])
    for id in paths.paths:
        total_l = 0
        for node in paths.getPathById(id):
#        print (paths.getPathById(id))
            if node in G.nodes:
                total_l += G.nodes[node]['length']
#short arm length condition
    #should we take in account already resolved rDNAs??
        if total_l < min_path_len or total_l > max_path_len:
            continue
        from_short_arm = False
        for node in paths.getPathById(id):
            for short_start in short_arm_starts:
                if short_start in dists and node in dists[short_start]:
                #path going from short arm telo to rdna to long arm
                    from_short_arm = True
        to_short_arm = False
        for node in rc_path(paths.getPathById(id)):
            for short_start in short_arm_starts:
                if node in dists and short_start in dists[node]:
                #path going from long arm to rdna to short arm
                    to_short_arm = True
        if from_short_arm and to_short_arm:
            print("Path {id} have orientation problems, skipping")
        else:
            if from_short_arm:
                res[id] = False
            elif to_short_arm:
                res[id] = True
    return res

#Map from path-id to direction (suspected telomere left, rdna right - true - otherwise false) 
def get_arm_paths_ids(telomere_locations_file, old_G, paths, rdna_nodes, min_path_len, max_path_len, ignore_telomeres):
    arm_paths_ids = {}
    aux_tel_nodes, G = get_telomeric_nodes(telomere_locations_file, old_G)

   # print(G.get_edge_data(">telomere_utig4-1871+_start"]))
    dists = dict(nx.all_pairs_dijkstra_path_length(G, weight='mid_length'))
    for id in paths.paths:
        total_l = 0
        for node in paths.getPathById(id):
#        print (paths.getPathById(id))
            if node in G.nodes:
                total_l += G.nodes[node]['length']
#short arm length condition
    #should we take in account already resolved rDNAs??
        if total_l < min_path_len or total_l > max_path_len:
            continue
        path = paths.getPathById(id)
        close_telo = False
        pathends = [path[0], path[-1]]
        start = 0
        telomeres_near = [False, False ]
        for pathend in pathends:
#TODO better use of coords            
#            if pathend in real_tel_nodes:
#                print (id + " in telomeres check")
#                close_telo = True
#oriented ones, so telnodes
            for tels in aux_tel_nodes:
                if start == 0:
                    order = [tels, pathend]
                else:
                    order = [pathend, tels]
                if order[0] in dists and order[1] in dists[order[0]] and dists[order[0]][order[1]] < G.nodes[pathend]['length'] + 200000:
                    print (f"{id} near telomeres check succedded, dir {start}")
                    close_telo = True
                    telomeres_near[start] = True    
            start +=1 
        close_rdna = False
        rdna_near = [False, False]
        cur_len = 0
        for node in path:
            if node in G.nodes:
                cur_len += G.nodes[node]['length']
            if node in rdna_nodes and cur_len < 200000:
                print (id + " rDNA in path start")
                close_rdna = True
                rdna_near[0] = True
                break
        path.reverse()
        cur_len = 0
        for node in path:
            if node in G.nodes:
                cur_len += G.nodes[node]['length']
            if node in rdna_nodes and cur_len < 200000:
                print (id + " rDNA in path end")
                close_rdna = True
                rdna_near[1] = True
                break
        path.reverse()

        direction = True
        closest_rdna_node = ""
        closest_rdna_dist = 999999999
        closest_suff_node = ""
        start = 0
        #need to check both directions...
        for suff in pathends:
            for t in rdna_nodes:
                if t == suff: 
                    continue
                #from rdna to paths' prefix, from paths' suffix to rdna 
                if start == 0: 
                    order = [t, suff]
                else:
                    order = [suff, t]
                if order[0] in dists and order[1] in dists[order[0]]:
                    if dists[order[0]][order[1]] < closest_rdna_dist:
                        closest_rdna_dist = dists[order[0]][order[1]] 
                        closest_rdna_node = t                        
                        closest_suff_node = suff
                    if dists[order[0]][order[1]] < G.nodes[suff]['length']+ G.nodes[t]['length'] + 2000000:
                        close_rdna = True
                        print (id + " rDNA near path")
                        rdna_near[start] = True 
                        break
            start += 1
        if not close_rdna and closest_rdna_dist < 999999999:
            print (id + " rDNA not so near path, closest node " + closest_rdna_node + " dist " + str(closest_rdna_dist) + " total_l " +str(G.nodes[closest_suff_node]['length']+ G.nodes[closest_rdna_node]['length']))
      
        if rdna_near[0] and rdna_near[1]:
            print (id + "  rDNA near both ends!!")
            if telomeres_near[0] == telomeres_near[1]:
                print (id + "  telomeres do not help, ignoring - no orientation")
                close_rdna = False
        else:
            if ignore_telomeres:
                #telomere can be near rdna side, we do not care about it for long arms
                telomeres_near[0] = True
                telomeres_near[1] = True

        if close_rdna and (close_telo or ignore_telomeres):
            if telomeres_near[0] and rdna_near[1]:
                print (id + " telomere->rdna")
                arm_paths_ids[id] = True
            elif telomeres_near[1] and rdna_near[0]:
                print (id + "rdna->telomere")
                arm_paths_ids[id] = False
            else:
                print (f"telo and rdna same side for {id}, ignoring")
                print (f"{id} telomeres {telomeres_near} rDNA {rdna_near} {ignore_telomeres}")
    return arm_paths_ids


def score(node, path, multiplicities, weights_map, use_multiple):
    res = 0
    for edge_ext in path:
        edge = edge_ext
        if edge in multiplicities and (multiplicities[edge] == 1 or use_multiple):
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
    if len(avg_links) == 0:
        return 0
    avg_links.sort()
    print (f"path {short_path} avg links {len(avg_links)} {avg_links[int(len(avg_links)/2)]}")
    return avg_links[int(len(avg_links)/2)]

#for components we allow to use nodes present in multiple paths, for haplo paths we do not
def scorepath(short_arm_path, long_path, mults, weights_map, use_multiple):    
    res = 0
    for node in short_arm_path:
        if node in mults and mults[node] == 1:
            res += score(node, long_path, mults, weights_map, use_multiple)
    return res

def path_length(path, G):
    res = 0
    for edge in path:
        if edge in G:
            res += G.nodes[edge]['length']
    return res

def get_best_path(short_id, path_container, paths, longarm_to_component, multiplicities, weights_map, G):
    path_scores = {}
    best_arms = []
    for long_arm in path_container.keys():
        path_scores[long_arm] = scorepath(paths.getPathById(short_id), path_container[long_arm], multiplicities, weights_map, True)
        best_arms.append(long_arm)
    best_arms.sort(key=lambda x: path_scores[x], reverse=True)    
    background_score = get_background_score(paths.getPathById(short_id), weights_map, G)
    print (short_id + " best arms " + str(best_arms) + " scores " + str([path_scores[x] for x in best_arms]))

    if len(best_arms) == 0:
        return "Unclear"
    if len(best_arms) == 1:
        if path_scores[best_arms[0]] > max (min(background_score * path_length(path_container[best_arms[0]], G)*3, 100), 10):
            return best_arms[0]
        else:
            print (f"One available path but score small or too close to background noise {path_scores[best_arms[0]]} and {background_score * path_length(path_container[best_arms[0]], G)}, skipping")
            return "Unclear"
    print (f"id {short_id}: best options are {best_arms[0]} and {best_arms[1]}, scores {path_scores[best_arms[0]]}  {path_scores[best_arms[1]]} background {background_score} {path_length(path_container[best_arms[0]], G) * background_score} {path_length(path_container[best_arms[1]], G) * background_score}")
    print (f"best updated scores {path_scores[best_arms[0]] - (path_length(path_container[best_arms[0]], G) * background_score)}  {(path_scores[best_arms[1]] - path_length(path_container[best_arms[1]], G) * background_score)}")
    print (f"best lengths  {path_length(path_container[best_arms[0]], G)}    {path_length(path_container[best_arms[1]], G)}")
    print (f"noise levels {(path_length(path_container[best_arms[0]], G) * background_score)} {(path_length(path_container[best_arms[1]], G) * background_score)}")
    
    if path_scores[best_arms[0]] <= max(min(background_score * path_length(path_container[best_arms[0]], G)*3, 100), 10):
        print (f"Best score is small or too close to background noise: {path_scores[best_arms[0]]} and {background_score * path_length(path_container[best_arms[0]], G)}, skipping")
        return "Unclear"
    if path_scores[best_arms[0]] - (path_length(path_container[best_arms[0]], G) * background_score) < 2 * (path_scores[best_arms[1]] - path_length(path_container[best_arms[1]], G) * background_score):
        print (f"No evident majority, best options are {best_arms[0]} and {best_arms[1]}, checking further")
        if longarm_to_component[best_arms[0]] != longarm_to_component[best_arms[1]]:
            #not homologous chromosomes, we can't resolve it
            return "Unclear"
        het_path_scores = {}
        for long_arm in path_container.keys():
            het_path_scores[long_arm] = scorepath(paths.getPathById(short_id), path_container[long_arm], multiplicities, weights_map, False)
        if het_path_scores[best_arms[0]] - (path_length(path_container[best_arms[0]], G) * background_score) < 2 * (het_path_scores[best_arms[1]] - path_length(path_container[best_arms[1]], G) * background_score):
            print(f"still no evident majority, best options are {best_arms[0]} and {best_arms[1]}, scores {path_scores[best_arms[0]]} and {path_scores[best_arms[1]]}, het scores {het_path_scores[best_arms[0]]} and {het_path_scores[best_arms[1]]}")
            return "Unclear"
    return best_arms[0]
