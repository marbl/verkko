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

# {node: path_id}
def get_nodes_in_unscaffed(all_paths, to_scaff):
    res = {}
    for id in to_scaff.keys():
        for node in all_paths.getPathById(id):
#            if node in multiplicities and multiplicities[node] == 1 and G.nodes[node]['length'] > long_enough:
            res[node] = id
    return res


# FIRST uniques where telomere is missing, {node: path_id}

def get_node_to_scaff(path, direction, multiplicities, G, long_enough=100000):
    if direction == "+":
        new_p = path
    else:
        new_p = reversed(path)
    for node in path:
        if node in multiplicities and multiplicities[node] == 1 and G.nodes[node]['length'] > long_enough:
            return node
    return "NONE"

def get_nodes_to_scaff(all_paths, to_scaff, multiplicities, G, long_enough=100000):
    res = {}
    for id in to_scaff.keys():
        path = all_paths.getPathById(id)
       
        if not (to_scaff[id][0]):
            node = get_node_to_scaff(path, "+", multiplicities, G, long_enough)
            for node in path:
                if node in multiplicities and multiplicities[node] == 1 and G.nodes[node]['length'] > long_enough:
                    res[node] = id
                    break
        if not (to_scaff[id][1]):    
            for node in reversed(path):
                if node in multiplicities and multiplicities[node] == 1 and G.nodes[node]['length'] > long_enough:
                    res[node] = id
                    break
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


def try_to_scaff(rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, indirectG, uncompressed_fasta):
    multiplicities = get_multiplicities(rukki_paths)
    tel_nodes, upd_G = get_telomeric_nodes(telomere_locations_file, G)
    to_scaff = get_paths_to_scaff(rukki_paths, tel_nodes, upd_G)
    nodes_unscaffed = get_nodes_in_unscaffed(rukki_paths, to_scaff)
    nodes_to_scaff = get_nodes_to_scaff(rukki_paths, to_scaff, multiplicities, G, 300000)

    coorded_connections = get_connections(hic_alignment_file, nodes_unscaffed)
    HiCGraph = graph_functions.loadHiCGraph(hic_alignment_file)
    matchGraph = graph_functions.loadMatchGraph(matches_file, indirectG, -239239239, 500000, 100000)
    lens = get_lengths(uncompressed_fasta)
    for node in nodes_to_scaff.keys():
        print (f"\nChecking extensions from {node}")
        connections = []
        shortened_node = node[:-1]
        if not shortened_node in HiCGraph.nodes():
            print (f"Node {node} not in hicGraph")
            continue
        if not node in nodes_unscaffed:
            print (f"Node {node} not in nodes_unscaffed")
            continue
        for nbr, datadict in HiCGraph.adj[shortened_node].items():
            connections.append([nbr, datadict['weight']])
        connections = sorted(connections, key=lambda x: x[1],reverse=True)

        for conn in connections:
            #no orientation from HiCGraph
            #for orientation in ("+", "-"):
            next_node = conn[0]
            next_path_id = "NONE"
            if next_node+"+" in nodes_unscaffed:
                next_path_id = nodes_unscaffed[next_node+"+"]
            if next_node+"-" in nodes_unscaffed:
                next_path_id = nodes_unscaffed[next_node+"-"]
            if [shortened_node, next_node] in matchGraph.edges and matchGraph.edges[shortened_node, next_node]['weight'] < 0:
                print (f"Skipping valid homology to {next_node}")
            elif next_path_id == nodes_unscaffed[node]:
                print (f"Skipping node {next_node} from same path")
            elif next_path_id != "NONE" and next_path_id != nodes_unscaffed[node]:
                print (f"Connection looks valid, from {node} to {next_node} {nodes_to_scaff[node]} {next_path_id}")
                print(get_pair_orientation((node[:-1], next_node), coorded_connections , lens))
                print (f"Path {nodes_unscaffed[node]} {rukki_paths.getPathString(nodes_unscaffed[node])}")
                print (f"aand {next_path_id} {rukki_paths.getPathString(next_path_id)}")
                if next_node+"+" in nodes_to_scaff or next_node+"-" in nodes_to_scaff:
                    print ("Two borderline top connected")
                else:
                    print ("NOT borderline top connected")
                break
            else:
                print (f"UNEXPECTED top connections from {node} to {next_node}")
                break

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
    
    CLOSE_TO_BORDER = 500000
    #we already know that htat node is first or last reliable enough node in the path
    def isPotentialNextPath(self, or_path_id, next_node):
        nor_path_id = or_path_id[:-1]
        path_orientation = or_path_id[-1]
        to_border = [0,0]
        path = self.rukki_paths.getPathById(nor_path_id)        
        for node in path:
            if node == next_node:
                break
            if not (node in self.compressed_lens):
                continue
            to_border[0] += self.compressed_lens[node]
        for node in reversed(path):
            if node == next_node:
                break
            if not (node in self.compressed_lens):
                continue
            to_border[1] += self.compressed_lens[node]
        sys.stderr.write (f"Checking path {or_path_id}, conn_node {next_node} dists to border {to_border} tels {self.scaffold_graph.nodes[or_path_id]['telomere']}\n")
        if path_orientation == "-":
            to_border.reverse()            
        
        if self.scaffold_graph.nodes[or_path_id]["telomere"][0]:
            sys.stderr.write (f" Trying to connect but telomere at the end!\n")            
            return False
        elif ((to_border[0] <= to_border[1]) or to_border[0] <= ScaffoldGraph.CLOSE_TO_BORDER):
            sys.stderr.write(f" Looks valid\n")
            return True
        else:
            sys.stderr.write(f"tried to connect, wrong direction\n")
            return False
    
    def fill_outgoing(self, nodes_unscaffed, or_path_id):
        #outgoing checks!
        print (or_path_id)
        print (self.scaffold_graph.nodes[or_path_id])
        if not (self.scaffold_graph.nodes[or_path_id]["telomere"][1]):
            path = self.rukki_paths.getPathById(or_path_id[:-1])
            direction = or_path_id[-1]
            or_last_node = get_node_to_scaff(path, direction, self.multiplicities, self.G, long_enough=500000)
            if or_last_node == "NONE":
                or_last_node = get_node_to_scaff(path, direction, self.multiplicities, self.G, long_enough=100000)
            if or_last_node == "NONE":
                sys.stderr.write(f"PATH {or_path_id} failed to find last node\n")
                return        
            if direction == "+":
                expected_orientation = or_last_node[-1]
            else:
                expected_orientation = self.rc_orientation(or_last_node[-1])
            nor_last_node = or_last_node[:-1]
            print (f"\nChecking extensions from {or_last_node}, expected best orientation {expected_orientation}")
            if not nor_last_node in self.HiCGraph.nodes():
                print (f"Node {nor_last_node} not in hicGraph")
                return
            if not or_last_node in nodes_unscaffed:
                print (f"Node {or_last_node} not in nodes_unscaffed")
                return
            connections = []                
            for nbr, datadict in self.HiCGraph.adj[nor_last_node].items():
                connections.append([nbr, datadict['weight']])
            connections = sorted(connections, key=lambda x: x[1],reverse=True)
            for conn in connections:
                #no orientation from HiCGraph
                #for orientation in ("+", "-"):
                next_node = conn[0]
                next_path_id = "NONE"
                #only one of those present, cause multiplicity= 1
                for next_orientation in ('-', '+'):
                    or_next_node = next_node + next_orientation
                    if or_next_node in nodes_unscaffed:
                        next_path_id = nodes_unscaffed[or_next_node]
                        break
                next_orientation = or_next_node[:-1]

                if [nor_last_node, next_node] in self.matchGraph.edges and self.matchGraph.edges[nor_last_node, next_node]['weight'] < 0:
                    print (f"Skipping valid homology to {next_node}")
                elif next_path_id == nodes_unscaffed[or_last_node]:
                    print (f"Skipping node {next_node} from same path")
                elif next_path_id != "NONE" and next_path_id != nodes_unscaffed[or_last_node]:
                    print (f"Connection looks valid, from {or_last_node} to {next_node} {or_path_id} {next_path_id}")
                    orientation_coints = get_pair_orientation((nor_last_node, next_node), self.coorded_connections , self.uncompressed_lens)
                    print (orientation_coints)
                    max_orient = max(orientation_coints, key=orientation_coints.get)
                    next_path_orientation = "+"
                    if max_orient[1] != next_orientation:
                        next_path_orientation = "-"
                    or_next_path_id = next_path_id + next_path_orientation
                    if self.isPotentialNextPath(or_next_path_id, or_next_node):
                        print (f"Path orientation check passed")
                        self.scaffold_graph.add_edge(or_path_id, or_next_path_id)
                    else:
                        print (f"Path orientation check failed!!")

#                    print (f"Path {nodes_unscaffed[last_or_node]} {rukki_paths.getPathString(nodes_unscaffed[last_or_node])}")
#                    print (f"aand {next_path_id} {rukki_paths.getPathString(next_path_id)}")
                    if next_node+"+" in self.nodes_to_scaff or next_node+"-" in self.nodes_to_scaff:
                        print ("Two borderline top connected")
                    else:
                        print ("NOT borderline top connected")
                    break
                else:
                    print (f"UNEXPECTED top connections from {nor_last_node} to {next_node}")
                    break
    
    def __init__(self, rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, uncompressed_fasta):
        self.multiplicities = get_multiplicities(rukki_paths)
        self.tel_nodes, self.upd_G = get_telomeric_nodes(telomere_locations_file, G)

        to_scaff = get_paths_to_scaff(rukki_paths, self.tel_nodes, self.upd_G)   
        nodes_unscaffed = get_nodes_in_unscaffed(rukki_paths, to_scaff)

#to be avoided
        self.nodes_to_scaff = get_nodes_to_scaff(rukki_paths, to_scaff, self.multiplicities, G, 300000)        
        self.coorded_connections = get_connections(hic_alignment_file, nodes_unscaffed)
 
        self.HiCGraph = graph_functions.loadHiCGraph(hic_alignment_file)
        self.matchGraph = graph_functions.loadMatchGraph(matches_file, G, -239239239, 500000, 100000)
        self.uncompressed_lens = get_lengths(uncompressed_fasta)
        self.compressed_lens = {}
        for node in G.nodes:
            self.compressed_lens[node] = G.nodes[node]['length']      
 
        self.rukki_paths = rukki_paths
        self.G = G

        self.scaffold_graph = nx.DiGraph()

        for id in rukki_paths.getPathIds():
            #connection between rc nodes?
            for dir in ('-', '+'):
                or_id = id + dir
                #TODO: what about shorter without  telomere
                tels_ends = [True, True]
                if id in to_scaff.keys():
                    if dir == '+':
                        tels_ends = to_scaff[id]
                    else:
                        tels_ends = [to_scaff[id][1], to_scaff[id][0]]
                self.scaffold_graph.add_node(or_id, telomere = tels_ends)    


        #possilby unefficient but whelp
        scores = {}
        all_connections = get_connections(hic_alignment_file, G.nodes)
        #naming, it's not outgoing!
        for outgoing_path_id in rukki_paths.getPathIds():
            scores[outgoing_path_id] = {}
            for incoming_path_id in rukki_paths.getPathIds():
                if incoming_path_id == outgoing_path_id:
                    continue
                #outgoing_path_id = "haplotype1_from_utig4-130"
                #incoming_path_id = "haplotype1_from_utig4-412"
                path_pair = [rukki_paths.getPathById(outgoing_path_id), rukki_paths.getPathById(incoming_path_id)]                
                scores[outgoing_path_id][incoming_path_id] = get_paths_connections(path_pair, all_connections, self.uncompressed_lens)
               
        for outgoing_path_id in rukki_paths.getPathIds():
            for incoming_path_id in rukki_paths.getPathIds():       
                if incoming_path_id == outgoing_path_id:
                    continue 
                if isHomologous(self.matchGraph, rukki_paths, [incoming_path_id, outgoing_path_id]):
                    continue
                for out_dir in ('-', '+'):
                    for in_dir in ('-', '+'):
                        or_out_id = outgoing_path_id + out_dir
                        or_in_id = incoming_path_id + in_dir

                        self.scaffold_graph.add_edge(or_out_id, or_in_id, weight = scores[outgoing_path_id][incoming_path_id][out_dir + in_dir])
                if outgoing_path_id in to_scaff and incoming_path_id in to_scaff:
                    print (f"Counted scores {outgoing_path_id} {incoming_path_id} {scores[outgoing_path_id][incoming_path_id]}")
        for outgoing_path_id in to_scaff:
            for out_dir in ('-', '+'):
                or_out_dir = outgoing_path_id + out_dir
                #telomere on the end - not interested
                if self.scaffold_graph.nodes[or_out_dir]['telomere'][1] == True:
                    continue

                local_scores = []
                for next_node in self.scaffold_graph.successors(or_out_dir):
                    local_scores.append([next_node, self.scaffold_graph.edges[or_out_dir, next_node]['weight']])
                print (f"Counting next from {or_out_dir}")
                local_scores.sort(key=lambda x: x[1], reverse=True)
                if len(local_scores) == 0:
                    print (f"Nothing next")  
                elif local_scores[0][1] <= 10:
                    print (f"very few links, best candidate {local_scores[0]}")                 
                elif len(local_scores) == 1:
                    print (f"Only next one, {local_scores[0]}")                
                elif local_scores[0][1] <  local_scores[1][1] * 1.5:
                    print (f"Not found next, first {local_scores[0]}, second best {local_scores[1]}")
                elif self.scaffold_graph.nodes[local_scores[0][0]]['telomere'][0]:
                    print (f"Best {local_scores[0]}, is good count but actually are in telomere!")
                else:
                    print (f"Really best {local_scores[0]}, second best {local_scores[1]}")

                

        
                #exit(0)
        '''
        for id in to_scaff:
            #connection between rc nodes?
            for dir in ('-', '+'):
                or_id = id + dir
                if dir == '+':
                    tels_ends = to_scaff[id]
                else:
                    tels_ends = [to_scaff[id][1], to_scaff[id][0]]
                
                self.scaffold_graph.add_node(or_id, telomere = tels_ends)    



        for or_path_id in self.scaffold_graph.nodes():
            print (f"checking node {or_path_id} {self.scaffold_graph.nodes[or_path_id]['telomere']}")
            self.fill_outgoing(nodes_unscaffed, or_path_id)
        '''

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

#TODO: this should actually be paths pair and not node pair!
#return scores for each of the orientations, ++, -+, +-, --,

#orientation within the path, it is not changing!
def get_pair_orientation(pair, orientations, connections, shift_before, shift_after, lens, middle):    
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
def get_paths_connections(paths, connections, lens, middle=2000000, ignore_short = 50000):
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
            scores_pair = get_pair_orientation(pair, orientations, connections, shift_before, shift_after, lens, middle)
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
