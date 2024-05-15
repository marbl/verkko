#!/usr/bin/env python3


import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions


#All functions expect directed nodes in utig4-234+ format, paths may contain gapped node

#path nodes stored with orientation, "utig4-234+" tsv format
class pathStorage:
    def __init__(self):        
        #from IDs to array of oriented+- nodes
        self.paths = {}
        #ignoring overlaps but who cares
        self.path_lengths = {}

    def getPathById(self, path_id):
        return self.paths[path_id]
    
    def getLength(self, path_id):
        return self.path_lengths[path_id]

    def getPathIds(self):
        return self.paths.keys()

#should not be used
    def getEdgeSequenceById(self, path_id):
        res = []
        for edge in self.paths[path_id]:
            #removing +- and gaps
            if edge[0] != "N":
                res.append(edge[:-1])
        return res
    def getPathString(self, path_id):
        return ",".join(self.paths[path_id])

    def addPath(self, line, G):
#        print (line)
        arr = line.strip().split()
        if len(arr) < 3:
            print (f"wrong path {line}")
            exit()
        separators = ">|<|,"
        edges = re.split(separators, arr[1])
        total_l = 0
        for edge in edges:
            node = edge
            if node in G.nodes:
                total_l += G.nodes[node]['length']
        edges = list(filter(None, edges))
#        print (edges)
        self.paths[arr[0]] = edges
        self.path_lengths[arr[0]] = total_l

    def addPathWithId(self, id, path, G):
        total_l = 0
        for edge in path:
            node = edge
            if node in G.nodes:
                total_l += G.nodes[node]['length']
        self.paths[id] = path
        self.path_lengths[id] = total_l
    
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
    path.reverse()
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
            for direction in ["+", "-"]:
                rdna_nodes.add(e + direction)
    return rdna_nodes

def get_telomeric_nodes(telomere_locations_file, G):
    aux_tel_nodes = set()
    new_G = G.copy()
    with open(telomere_locations_file) as f:
        for l in f:
            parts = l.strip().split('\t')
            telnode = ""
            graph_node = parts[0]
            if int(parts[1]) < 20000:
                telnode="telomere_"+graph_node+"+_start"
                new_G.add_node(telnode, length=0, coverage=0)
                new_G.add_edge(telnode, graph_node+'+', mid_length=G.nodes[graph_node+'+']['length'])
                aux_tel_nodes.add(telnode)

                telnode="telomere_"+graph_node+"-_end"
                new_G.add_node(telnode, length=0, coverage=0)
                new_G.add_edge(graph_node+'-', telnode, mid_length=G.nodes[graph_node+'+']['length'])
                aux_tel_nodes.add(telnode)

            if int(parts[2])+20000 > G.nodes[graph_node + '+']['length']:
                telnode=">telomere_"+graph_node+"+_end"
                new_G.add_node(telnode, length=0, coverage=0)
                new_G.add_edge(graph_node+'+', telnode, mid_length=G.nodes[graph_node+'+']['length'])  
                aux_tel_nodes.add(telnode)

                telnode="telomere_"+graph_node+"-_start"
                new_G.add_node(telnode, length=0, coverage=0)
                new_G.add_edge(telnode, graph_node+'-', mid_length=G.nodes[graph_node+'+']['length'])     
                aux_tel_nodes.add(telnode)

    return aux_tel_nodes, new_G

def read_rukki_paths(rukki_file, G):
    res = pathStorage()
    for line in open(rukki_file):
        arr = line.strip().split()
        if arr[0] == "name":
            continue
        res.addPath(line.strip(), G)
    return res

def get_multiplicities(paths):
    multiplicities = {}
    for path_id in paths.getPathIds():
        for edge in paths.getEdgeSequenceById(path_id):
            for dir in ["+", "-"]:
                node = edge + dir      
                if not node in multiplicities:
                    multiplicities[node] = 0
                multiplicities[node] += 1
    return multiplicities

#returns dict, {id:[present_start_relo, present_end_telo]}
def get_paths_to_scaff(paths, tel_nodes, dirG, long_enough=500000):
    res = {}
    for id in paths.paths:
        total_l = paths.path_lengths[id]
        tel_start = False
        tel_end = False    
        for telomere in tel_nodes:
            if dirG.has_edge(telomere, paths.paths[id][0]):
                tel_start = True
            if dirG.has_edge(paths.paths[id][-1], telomere):
                tel_end = True
        if tel_start and tel_end:
            continue
        #all long enough AND containing telomere
        if total_l > long_enough:
            res[id] = [tel_start, tel_end]
            print (f"will use path {paths.paths[id]} {tel_start} {tel_end}")
    return res

def get_lengths(fasta_file):
    res = {}
    cur = ""
    for line in open(fasta_file):
        if line[0] == ">":
            cur = line.strip()[1:]
        else:
            res[cur] = len(line.strip())
            cur = ""
    return res

def isHomologous (match_graph, paths, ids):
    hom_size = 0
    for p0 in paths.getPathById(ids[0]):
        for p1 in paths.getPathById(ids[1]):
            nor_p0 = p0.strip("-+")
            nor_p1 = p1.strip("-+")
            if match_graph.has_edge(nor_p0, nor_p1):
                if match_graph.edges[nor_p0, nor_p1]['weight'] < 0:
                    
                    hom_size += match_graph.edges[nor_p0, nor_p1]['homology_len']
    if hom_size * 2> paths.getLength(ids[0])  or hom_size * 2> paths.getLength(ids[1]):
        print (f"Found homo paths {ids}")
        return True
    else:
        return False
#TODO: or inherit from nx.Digraph??
class ScaffoldGraph:

    def rc_orientation(self, c):
        if c == "+":
            return "-"
        if c == "-":
            return "+"
        return c    

    def rc_path_id(self, path_id):
        return path_id[:-1] + self.rc_orientation(path_id[-1])
    
    def getHaploidPaths(self):
        haploids = set()
        nodes_to_path_ids = {}        
        for nor_path_id in self.rukki_paths.getPathIds():
            for or_node in self.rukki_paths.getPathById(nor_path_id):
                nor_node = or_node.strip("-+")
                if not (nor_node in nodes_to_path_ids):
                    nodes_to_path_ids[nor_node] = []
                nodes_to_path_ids[nor_node].append(nor_path_id)
        #possibly we'll need that graph but not now
        #homGraph = nx.Graph()
        for nor_path_id in self.rukki_paths.getPathIds():
            #let's leave for the graph.
            homs = {}
            total_hom = 0
            for or_node in self.rukki_paths.getPathById(nor_path_id):
                nor_node = or_node.strip("-+")
                if not (nor_node in self.matchGraph.nodes):
                    continue
                for next in self.matchGraph.neighbors(nor_node):
                    if self.matchGraph.edges[nor_node, next]['weight'] < 0:
                        if len (nodes_to_path_ids[next]) > 1:
                            #soemthing weird, ignoring
                            continue
                        next_id = nodes_to_path_ids[next][0]
                        if not (next_id in homs):
                            homs[next_id] = 0
                        homs[next_id] += self.matchGraph.edges[nor_node, next]['homology_len']
                        total_hom += self.matchGraph.edges[nor_node, next]['homology_len']
            path_len = self.rukki_paths.getLength(nor_path_id)
            if total_hom * 2 < path_len:
                haploids.add(nor_path_id)
                if path_len > 2000000:
                    print (f"Found haploid path {nor_path_id} with homology {total_hom} and len {path_len} ")
        return haploids
    
    def __init__(self, rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, uncompressed_fasta):
        self.multiplicities = get_multiplicities(rukki_paths)
        self.tel_nodes, self.upd_G = get_telomeric_nodes(telomere_locations_file, G)
        self.to_scaff = get_paths_to_scaff(rukki_paths, self.tel_nodes, self.upd_G)   
        self.hic_alignment_file = hic_alignment_file
        self.matchGraph = graph_functions.loadMatchGraph(matches_file, G, -239239239, 500000, 100000)
        self.uncompressed_lens = get_lengths(uncompressed_fasta)
        self.compressed_lens = {}
        for node in G.nodes:
            self.compressed_lens[node] = G.nodes[node]['length']      
            self.compressed_lens[node.strip("-+")] = G.nodes[node]['length']      

        self.rukki_paths = rukki_paths
        self.G = G
        self.dists = dict(nx.all_pairs_dijkstra_path_length(G, weight=lambda u, v, d: G.edges[u, v]['mid_length']))

        self.haploids = self.getHaploidPaths()
        self.scaffold_graph = nx.DiGraph()

        for id in rukki_paths.getPathIds():
            for dir in ('-', '+'):
                or_id = id + dir
                #TODO: what about shorter without  telomere
                tels_ends = [True, True]
                if id in self.to_scaff.keys():
                    if dir == '+':
                        tels_ends = self.to_scaff[id]
                    else:
                        tels_ends = [self.to_scaff[id][1], self.to_scaff[id][0]]
                self.scaffold_graph.add_node(or_id, telomere = tels_ends)    
        #possilby unefficient but whelp
        scores = {}
        all_connections = self.get_connections(self.hic_alignment_file)
        for from_path_id in self.rukki_paths.getPathIds():
            scores[from_path_id] = {}
            for to_path_id in self.rukki_paths.getPathIds():
                if to_path_id == from_path_id:
                    continue               
                scores[from_path_id][to_path_id] = self.getPathPairConnections([from_path_id, to_path_id], all_connections, self.uncompressed_lens)
                
        for from_path_id in rukki_paths.getPathIds():
            for to_path_id in rukki_paths.getPathIds():       
                if to_path_id == from_path_id:
                    continue 
                if isHomologous(self.matchGraph, rukki_paths, [to_path_id, from_path_id]):
                    continue
                for from_dir in ('-', '+'):
                    for to_dir in ('-', '+'):
                        or_from_path_id = from_path_id + from_dir
                        or_to_path_id = to_path_id + to_dir
                        self.scaffold_graph.add_edge(or_from_path_id, or_to_path_id, weight = scores[from_path_id][to_path_id][from_dir + to_dir])
                if from_path_id in self.to_scaff and to_path_id in self.to_scaff:
                    print (f"Counted scores {from_path_id} {to_path_id} {scores[from_path_id][to_path_id]}")

    def forbiddenPair(self, from_path_id, to_path_id):    
        nor_from_path_id = from_path_id.strip('-+')
        nor_to_path_id = to_path_id.strip('-+')
        #Heterogametic chromosomes get more links since there is no homologous one to absorb multiple alignments, so no connection of diploid and long enough ahploids    
        if nor_from_path_id in self.haploids and self.rukki_paths.getLength(nor_from_path_id) > 3000000 and not (nor_to_path_id in self.haploids):
            return True
        if nor_to_path_id in self.haploids and self.rukki_paths.getLength(nor_to_path_id) > 3000000 and not (nor_from_path_id in self.haploids):
            return True
        #relatively short fragments with telomere are special case, we may fail to detect orientation there but belive in telomere.
        if self.rukki_paths.getLength(nor_to_path_id) <= 3000000 and self.scaffold_graph.nodes[to_path_id]['telomere'][0]:
            return True
        
        return False
    
    #Main logic is here!        
    def findExtension(self, cur_path_id):
        local_scores = []
        print (f"Checking {cur_path_id}")
        if self.scaffold_graph.nodes[cur_path_id]['telomere'][1]:
            print (f"Stopped at the telomere")
            return "NONE"
        for next_node in self.scaffold_graph.successors(cur_path_id):
            #specific hacks to avoid
            local_scores.append([next_node, self.scaffold_graph.edges[cur_path_id, next_node]['weight']])            

        local_scores.sort(key=lambda x: x[1], reverse=True)
        best_ind = -1
        second_best_ind = -1        
        for i in range (0, len(local_scores)):
            if not self.forbiddenPair(cur_path_id, local_scores[i][0]):
                best_ind = i
                break
        for j in range (i+1, len(local_scores)):
            if not self.forbiddenPair(cur_path_id, local_scores[j][0]):
                second_best_ind = j
                break
        if len(local_scores) == 0:            
            print (f"Nothing next")  
            return "NONE"
        elif local_scores[best_ind][1] <= 10:
            print (f"very few links, best valid candidate {local_scores[best_ind]}")                 
            return "NONE"
        #not valid but waay best solution exists
        elif len(local_scores) == 1:
            print (f"Only next one, {local_scores[best_ind]}") 
            return local_scores[best_ind][0]                       
        elif local_scores[best_ind][1] <  local_scores[second_best_ind][1] * 1.5:
            print (f"Not found next, first {local_scores[best_ind]}, second best {local_scores[second_best_ind]}")
            return "NONE"
        elif self.scaffold_graph.nodes[local_scores[best_ind][0]]['telomere'][0]:
            print (f"Best {local_scores[best_ind]}, is good count but actually are in telomere!")            
            return "NONE"
        else:
            print (f"Really best {local_scores[best_ind]}, second best {local_scores[second_best_ind]}")            
            return local_scores[best_ind][0]
        
    def generateScaffolds(self):
        res = []
        #will grow to right these paths in length order 
        tel_starting_paths = []        
        middle_paths = []
        #to avoid outputing same path twice
        nor_used_path_ids = set()
        for from_path_id in self.to_scaff:
            for from_dir in ('-', '+'):
                or_from_path_id = from_path_id + from_dir
                if not self.scaffold_graph.nodes[or_from_path_id]['telomere'][1]:
                    if self.scaffold_graph.nodes[or_from_path_id]['telomere'][0]:
                        tel_starting_paths.append(or_from_path_id)
                    else:
                        middle_paths.append(or_from_path_id)        

        tel_starting_paths.sort(key=lambda x: self.rukki_paths.getLength(x.strip('-+')), reverse=True)        
        middle_paths.sort(key=lambda x: self.rukki_paths.getLength(x.strip('-+')), reverse=True)
        starting_paths = tel_starting_paths + middle_paths
        print ("Starting paths")
        print (starting_paths)
        for or_from_path_id in starting_paths:
            if or_from_path_id.strip('-+') in nor_used_path_ids:
                continue            
            cur_scaffold = [or_from_path_id]
            cur_path_id = or_from_path_id            
            nor_used_path_ids.add(or_from_path_id.strip('-+'))
            while True:
                next_path_id = self.findExtension(cur_path_id)
                if next_path_id == "NONE":
                    break
                elif next_path_id.strip('-+') in nor_used_path_ids:
                    print (f"Extention {next_path_id} looks good but already used")
                    break
                if self.rc_path_id(self.findExtension(self.rc_path_id(next_path_id))) != cur_path_id:
                    print (f"backward check failed for {next_path_id}")
                    break
                print (f"Extending {cur_path_id} with {next_path_id}")
                cur_scaffold.append(next_path_id)
                nor_used_path_ids.add(next_path_id.strip('-+'))
                cur_path_id = next_path_id
            print (f"scaffold {cur_scaffold} \n")
            res.append(cur_scaffold)
        return res

    #returns: dict {(start_id, end_id):[[start_pos1, end_pos1]]}. Coords not compressed!
    def get_connections(self, alignment_file):
        res = {}
        #A01660:39:HNYM7DSX3:1:1101:1696:29982   utig4-73        utig4-1056      1       16949880        78591191
        for line in open (alignment_file):
            arr = line.split()
            if not (arr[1], arr[2]) in res:
                res[(arr[1], arr[2])] = []
#            print (arr[1], arr[2])
#            print (arr)
            next = [int(arr[4]), int(arr[5])]
            res[(arr[1], arr[2])].append(next)
            if not (arr[2], arr[1]) in res:
                res[(arr[2], arr[1])] = []
            res[(arr[2], arr[1])].append([int(arr[5]), int(arr[4])])  
        return res

    #return scores for each of the orientations, ++, -+, +-, --,
    #orientation within the path, it is not changing!
    def getNodePairConnections(self, pair, orientations, connections, shift_before, shift_after, lens, middle):    
        #This is scores for PATH orientation
        scores = {"++":0, "-+":0, "+-":0, "--":0, "middle":0}    
        #pair = ("utig4-1799", "utig4-1957")    
#        print (self.matchGraph.edges[pair[0], pair[1]]['homology_len'])
#        print (self.matchGraph.edges[pair[0], pair[1]]['intervals'])
#        exit (0)
        filtered = 0
        not_filtered = 0
        if self.matchGraph.has_edge(pair[0], pair[1]):
            intervals = self.matchGraph.edges[pair[0], pair[1]]['intervals']
        else:
            intervals = [[],[]]
        #print (intervals)
        for conn in connections[pair]:        
            #filter "middle" nodes
            in_homo = False
            for i in range (0, len(intervals[0])):
                local_homo = True
                for j in range (0, 2):
                    if conn[j] <  intervals[j][i][0] or conn[j] > intervals[j][i][1]:
                        local_homo = False
                if not local_homo:
                    in_homo = True
                    break
            if in_homo:
                filtered+= 1
                continue
            else:
                not_filtered += 1

            near_ends = [False, False]            
            dists_to_end =[[0,0],[0,0]]
            for i in range (0, 2):
                fixed_coords = [conn[i], lens[pair[i]] - conn[i]]
                if orientations[i] == "-":
                    fixed_coords.reverse()            
                dists_to_end[i][0] = shift_before[i]
                dists_to_end[i][1] = shift_after[i]
                for j in range(0, 2):
                    dists_to_end[i][j] += fixed_coords[j]
                for j in range (0, 2):
                    if (dists_to_end[i][j]< middle):
                        near_ends[i] = True        
            str = ""
            if dists_to_end[0][0] < dists_to_end[0][1]:
                str += "-"
            else:
                str += "+"
            if dists_to_end[1][0] > dists_to_end[1][1]:
                str += "-"
            else:
                str += "+"
            if near_ends[0] and near_ends[1]:
                scores[str] += 1
            else:
                scores["middle"] += 1
        print (f"scores for {pair} {scores} {orientations}")
        print (f"shifts {shift_before} {shift_after}")
        print (f" {intervals}")
        print (f" {pair} filtered/not_filtered {filtered} {not_filtered}")
        #exit(0)
        return scores
    
    #TODO: move to matchGraph
    def homologousOrNodes (self, or_node):
        nor_node = or_node.strip("-+")
        if not nor_node in self.matchGraph.nodes:
            return set()
        orient = or_node[-1]
        res = set()
        for hom_node in self.matchGraph.neighbors(nor_node):
            if self.matchGraph.edges[nor_node, hom_node]['weight'] < 0:
                if self.matchGraph.edges[nor_node, hom_node]['orientation'] == '+':
                    res.add(hom_node + orient)
                else:
                    res.add (hom_node+ self.rc_orientation(orient))
        return res
    
    def fixOrientation(self, path_ids, scores):
        #orientation for telomere-containing short contigs is fixed, but because of filtering (i.e. distal bits) read counts can be misleading
        #first element should _not_ swapped for telomere on the left, second should.
        correct_orientations="+-"
        #path_ids = ["haplotype2_unused_utig4-1963", "haplotype2_unused_utig4-465"]
        total_score = 0
        for i in ('-', '+'):
            for j in ('-', '+'):
                total_score += scores[i+j]
        if total_score > 0:
            for i in range (0, 2):
                #Do we need length check here at all?
                if self.rukki_paths.getLength(path_ids[i]) <= 30000000:
                    if self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1]:
                        correct_or = correct_orientations[i]
                    elif self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0]:
                        correct_or = self.rc_orientation(correct_orientations[i])
                    else:
                        #connectivity check

                        #alternative path node set, all orientations
                        alt_nodes = set()
                        for node in self.rukki_paths.getEdgeSequenceById(path_ids[1 - i]):
                            for orient in ('-', '+'):       
                                if not node+orient in self.G.nodes:
                                    continue                         
                                alt_nodes.add(node + orient)
                                #there can be gap in one of the haplotypes, so using oriented nodes  
                                for hom_node in self.homologousOrNodes(node + orient):
                                    alt_nodes.add(hom_node)
                        check_nodes = {'-':set(), '+':set()}
                        for node in self.rukki_paths.getPathById(path_ids[i]):
                            #current_node orientation
                            if not node in self.G.nodes:
                                continue
                            for node_orient in ('-', '+'): 
                            #different orientations of the PATH!                                    
                                if node_orient == node[-1]:                                    
                                    #no rc to path!
                                    path_orient = "+"
                                else:
                                    path_orient = "-"
                                to_check = node.strip("-+") + node_orient
                                check_nodes[path_orient].add(to_check)
                                #again, to work with gaps    
                                for hom_node in self.homologousOrNodes(to_check):
                                    check_nodes[path_orient].add(hom_node)

                        print (f" path_pair {path_ids} counting dists") 
                        shortest_paths = {'-':1000000000, '+':1000000000}
                        for orient in ('-', '+'):
                            for alt_node in alt_nodes:
                                for check_node in check_nodes[orient]:
                                    check_dist = 1000000000
                                    if i == 0 and alt_node in self.dists[check_node]:
                                        check_dist = self.dists[check_node][alt_node] 
                                    elif i == 1 and check_node in self.dists[alt_node]: 
                                        check_dist = self.dists[alt_node][check_node]
                                    if check_dist == 1000000000:
                                        continue
                                #distance is between centers of the nodes, G is compressed:(
                                    tuned_dist = (check_dist - self.compressed_lens[check_node.strip('-+')] - self.compressed_lens[alt_node.strip('-+')])/ 2
                                    #print (f"Pairlens: {alt_node}  {check_node} {i} {check_dist} {self.compressed_lens[check_node.strip('-+')]} {self.compressed_lens[alt_node.strip('-+')]} {tuned_dist}")

                                    check_dist = tuned_dist
                                    if check_dist < shortest_paths[orient]:
                                        shortest_paths[orient] = check_dist
                        print (f" path_pair {path_ids} swapping {i} shortest paths {shortest_paths}")
                        #One orientation is close in graph, second is not
                        min_cutoff = min(500000, self.rukki_paths.getLength(path_ids[i]) / 4)                        
                        max_cutoff = self.rukki_paths.getLength(path_ids[i]) * 3 / 4
                        if shortest_paths['-'] < min_cutoff and shortest_paths['+'] > max_cutoff:
                            correct_or = "-"    
                        elif shortest_paths['+'] < min_cutoff / 4 and shortest_paths['-'] > max_cutoff:
                            correct_or = "+"
                        else:
                            correct_or = ""
                    if correct_or != "":            
                        pair_or = ["",""]
                        pair_or[i] = correct_or
                        
                        print (f"tuning pair {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                        for second_or in ('-', '+'):
                            correct_pair = pair_or.copy()                       
                            correct_pair[1 - i] = second_or
                            incorrect_pair = correct_pair.copy()
                            incorrect_pair[i] = self.rc_orientation(correct_or)
                            correct_pair = "".join(correct_pair)
                            incorrect_pair = "".join(incorrect_pair)
                            print (f"moving {incorrect_pair} to {correct_pair}")
                            scores[correct_pair] += scores[incorrect_pair]
                            scores[incorrect_pair] = 0
                        print (f"tuned pair {path_ids}, scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
        
        return scores
    
    #TODO: this should actually be paths pair and not node pair!
    #return scores for each of the orientations, ++, -+, +-, --,
    def getPathPairConnections(self, path_ids, connections, lens, middle=5000000, ignore_short = 50000):
        #from start to occurance of node exclusive, do not care about multiplicity > 1 (possibly should filter)
        length_before = [{}, {}]
        total_len = [0, 0]
        paths = [self.rukki_paths.getPathById(path_ids[0]), self.rukki_paths.getPathById(path_ids[1])]
        for i in range (0, 2):
            for or_node in paths[i]:
                length_before[i][or_node.strip('-+')] = total_len[i]
                if or_node.strip('-+') in lens:
                    total_len[i] += lens[or_node.strip('-+')]
        scores = {"++":0, "-+":0, "+-":0, "--":0, "middle":0}
        for first in paths[0]:
            nor_f = first.strip('-+')
            if not (nor_f in lens) or lens[nor_f] < ignore_short:
                continue
            for second in paths[1]:
                nor_s = second.strip('-+')
                if not (nor_s in lens) or lens[nor_s] < ignore_short:
                    continue
                pair = (nor_f, nor_s)
                orientations = (first[-1], second[-1])
                if not pair in connections:
                    continue
                shift_before = []
                shift_after = []
                for i in range (0, 2):
                    shift_before.append(length_before[i][pair[i]])
                    shift_after.append(total_len[i] - shift_before[i] - lens[pair[i]])
                scores_pair = self.getNodePairConnections(pair, orientations, connections, shift_before, shift_after, lens, middle)
                for key in scores_pair:
                    scores[key] += scores_pair[key] 
        scores = self.fixOrientation(path_ids, scores)
        return scores

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
