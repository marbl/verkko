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
        self.rukki_paths = rukki_paths
        self.G = G
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
        all_connections = get_connections(self.hic_alignment_file, self.G.nodes)
        for from_path_id in self.rukki_paths.getPathIds():
            scores[from_path_id] = {}
            for to_path_id in self.rukki_paths.getPathIds():
                if to_path_id == from_path_id:
                    continue
                path_pair = [rukki_paths.getPathById(from_path_id), rukki_paths.getPathById(to_path_id)]                
                scores[from_path_id][to_path_id] = getPathPairConnections(path_pair, all_connections, self.uncompressed_lens)
               
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
        starting_paths = []        
        #to avoid outputing same path twice
        nor_used_path_ids = set()
        for from_path_id in self.to_scaff:
            for from_dir in ('-', '+'):
                or_from_path_id = from_path_id + from_dir
                if not self.scaffold_graph.nodes[or_from_path_id]['telomere'][1]:
                    starting_paths.append(or_from_path_id)

        starting_paths.sort(key=lambda x: self.rukki_paths.getLength(x.strip('-+')), reverse=True)        
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
            res.append(cur_scaffold)
        return res

#returns: dict {(start_id, end_id):[[start_pos1, end_pos1]]}. Coords not compressed!
def get_connections(alignment_file, interesting_nodes):
    res = {}
    undirected_interesting = {}
    for node in interesting_nodes:
        undirected_interesting[node[:-1]] = interesting_nodes[node]
    #A01660:39:HNYM7DSX3:1:1101:1696:29982   utig4-73        utig4-1056      1       16949880        78591191
    for line in open (alignment_file):
        arr = line.split()
        if arr[1] in undirected_interesting and arr[2] in undirected_interesting:
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
def getNodePairConnections(pair, orientations, connections, shift_before, shift_after, lens, middle):    
    scores = {"++":0, "-+":0, "+-":0, "--":0, "middle":0}
    for conn in connections[pair]:        
        #filter "middle" nodes
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
        if conn[0] < lens[pair[0]] - conn[0]:
            str += "-"
        else:
            str += "+"
        if conn[1] > lens[pair[1]] - conn[1]:
            str += "-"
        else:
            str += "+"
        if near_ends[0] and near_ends[1]:
            scores[str] += 1
        else:
            scores["middle"] += 1
    #print (f"scores for {pair} {scores} {orientations}")
    #print (f"shifts {shift_before} {shift_after}")
    return scores
                
#TODO: this should actually be paths pair and not node pair!
#return scores for each of the orientations, ++, -+, +-, --,
def getPathPairConnections(paths, connections, lens, middle=2000000, ignore_short = 50000):
    #from start to occurance of node exclusive, do not care about multiplicity > 1 (possibly should filter)
    length_before = [{}, {}]
    total_len = [0, 0]
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
            scores_pair = getNodePairConnections(pair, orientations, connections, shift_before, shift_after, lens, middle)
            for key in scores_pair:
                scores[key] += scores_pair[key] 

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
