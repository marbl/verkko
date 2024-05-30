#!/usr/bin/env python3

import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions
import logging
from scaffolding import logger_wrap, match_graph, path_storage
import rdna_scaff_functions as sf

#TODO: or inherit from nx.Digraph??
class ScaffoldGraph:
    LONG_HAPLOID_CUTOFF = 5000000
    #should not be actuallly used, since they'll be reordered with MAX_REORDERING_LENGTH check
    SHORT_TEL_CUTOFF = 5000000
    MIN_LINKS = 10
    CLEAR_MAJORITY = 1.5    
    MAX_REORDERING_LENGTH = 30000000

    #If not near end, connection is ignored
    NEAR_PATH_END = 5000000
    SHORT_INGORED_NODE = 20000

    CLOSE_IN_GRAPH = 500000

    MATCHGRAPH_LONG_NODE = 500000
    MATCHGRAPH_MIN_ALIGNMENT = 100000

    def __init__(self, rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, uncompressed_fasta, logger):
        self.logger = logger_wrap.UpdatedAdapter(logger, self.__class__.__name__)
        self.multiplicities = sf.get_multiplicities(rukki_paths)
        #upd_G - graph with telomeric nodes
        self.tel_nodes, self.upd_G = sf.get_telomeric_nodes(telomere_locations_file, G)
        self.logger.debug("Telomeric nodes")
        self.logger.debug(self.tel_nodes)
        self.output_basename = "scaff_rukki.paths"

        #Used for scaffolding starts and debug, 
        self.to_scaff = sf.get_paths_to_scaff(rukki_paths, self.tel_nodes, self.upd_G)   
        self.hic_alignment_file = hic_alignment_file

        #TODO: duplicated... self.matchGraph should not be refered directly
        self.mg = match_graph.MatchGraph(matches_file, G, -239239239, ScaffoldGraph.MATCHGRAPH_LONG_NODE, ScaffoldGraph.MATCHGRAPH_MIN_ALIGNMENT, logger)
        self.matchGraph = self.mg.getMatchGraph()
        
        #debug purposes
        self.dangerous_swaps = {}

        self.uncompressed_lens = sf.get_lengths(uncompressed_fasta)
        self.compressed_lens = {}
        for node in self.upd_G.nodes:
            self.compressed_lens[node] = self.upd_G.nodes[node]['length']      
            self.compressed_lens[node.strip("-+")] = self.upd_G.nodes[node]['length']  
            

        self.rukki_paths = rukki_paths
        self.G = G
        self.dists = dict(nx.all_pairs_dijkstra_path_length(self.upd_G, weight=lambda u, v, d: self.upd_G.edges[u, v]['mid_length']))
        self.logger.info("Pairwise distances in assembly graph calculated")
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
                if self.mg.isHomologous([rukki_paths.getPathById(from_path_id), rukki_paths.getPathById(to_path_id)], [rukki_paths.getLength(from_path_id), rukki_paths.getLength(to_path_id)]):
                    continue
                for from_dir in ('-', '+'):
                    for to_dir in ('-', '+'):
                        or_from_path_id = from_path_id + from_dir
                        or_to_path_id = to_path_id + to_dir
                        self.scaffold_graph.add_edge(or_from_path_id, or_to_path_id, weight = scores[from_path_id][to_path_id][from_dir + to_dir])
                if from_path_id in self.to_scaff and to_path_id in self.to_scaff:
                    self.logger.debug (f"Counted scores {from_path_id} {to_path_id} {scores[from_path_id][to_path_id]}")

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
                    self.logger.info(f"Found haploid path {nor_path_id} with homology {total_hom} and len {path_len} ")
        return haploids


    def getClosestTelomere(self, path, direction):
        #From telomere to rc(path_end)
        if direction == '+':
            path = sf.rc_path(path)
        closest = 1000000000

        add_dist = 0
        #Or use all nodes?
        shortened_path = [path[0]]

        for tel_node in self.tel_nodes:
            closest = min(closest, self.nodeToPathDist(tel_node, shortened_path, False) + add_dist)
        return closest/2
    
    #Dist from node end to path. Allowed to go not in the first path node, but then additional length in path is added
    #Optionally allowing to use homologous nodes (to improve in gaps)
    def nodeToPathDist(self, node, path, check_homologous):
        closest = 1000000000
        add_dist = 0
        for path_node in path:
            if path_node in self.upd_G.nodes:
                if path_node in self.dists[node]:
                    if path_node == node:
                        closest = 0
                    closest = min(closest, self.dists[node][path_node] + add_dist - self.compressed_lens[node.strip("-+")] - self.compressed_lens[path_node.strip("-+")])
                    if check_homologous:
                        for hom_node in self.homologousOrNodes(path_node):
                            if hom_node == node:
                                closest = 0
                            if hom_node in self.dists[node]:
                               closest = min(closest, self.dists[node][hom_node] + add_dist - self.compressed_lens[node.strip("-+")] - self.compressed_lens[hom_node.strip("-+")])
                add_dist += 2* self.compressed_lens[path_node.strip("-+")]
        return closest/2
    #Optionally allowing to use homologous nodes (to improve in gaps)

    def pathDist(self, path_from, path_to, check_homologous):
        closest = 1000000000
        add_dist = 0
        for node in reversed(path_from):
            if node in self.upd_G.nodes:
                closest = min(closest, self.nodeToPathDist(node, path_to, check_homologous) + add_dist)
                if check_homologous:
                    for hom_node in self.homologousOrNodes(node) :
                        closest = min(closest, (self.nodeToPathDist(hom_node, path_to, check_homologous) + add_dist))
                add_dist += 2* self.compressed_lens[node.strip("-+")]
        return closest

            
    def forbiddenPair(self, from_path_id, to_path_id):    
        nor_from_path_id = from_path_id.strip('-+')
        nor_to_path_id = to_path_id.strip('-+')
        #Heterogametic chromosomes get more links since there is no homologous one to absorb multiple alignments, so no connection of diploid and long enough ahploids    
        if nor_from_path_id in self.haploids and self.rukki_paths.getLength(nor_from_path_id) > ScaffoldGraph.LONG_HAPLOID_CUTOFF and not (nor_to_path_id in self.haploids):
            self.logger.warning(f"Banning link from long haploid {from_path_id} to diploid {to_path_id}")
            return True
        if nor_to_path_id in self.haploids and self.rukki_paths.getLength(nor_to_path_id) > ScaffoldGraph.LONG_HAPLOID_CUTOFF and not (nor_from_path_id in self.haploids):
            self.logger.warning(f"Banning link from diploid {from_path_id} to long haploid {to_path_id}")
            return True
        #relatively short fragments with telomere are special case, we may fail to detect orientation there but belive in telomere.
        if self.rukki_paths.getLength(nor_to_path_id) <= ScaffoldGraph.SHORT_TEL_CUTOFF and self.scaffold_graph.nodes[to_path_id]['telomere'][0]:
            self.logger.error(f"Banning link from {from_path_id} into short telomere, should not happen {to_path_id}")
            return True
        return False
    
    #Main logic is here!        
    def findExtension(self, cur_path_id):
        local_scores = []
        self.logger.info(f"Checking {cur_path_id}")
        if self.scaffold_graph.nodes[cur_path_id]['telomere'][1]:
            self.logger.info (f"Stopped at the telomere")
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
            self.logger.info (f"Nothing next")  
            return "NONE"
        elif local_scores[best_ind][1] <= ScaffoldGraph.MIN_LINKS:
            self.logger.info (f"very few links, best valid candidate {local_scores[best_ind]}")                 
            return "NONE"
        #not valid but waay best solution exists
        elif len(local_scores) == 1:
            self.logger.info (f"Only next one, {local_scores[best_ind]}") 
            return local_scores[best_ind][0]                       
        elif local_scores[best_ind][1] <  local_scores[second_best_ind][1] * ScaffoldGraph.CLEAR_MAJORITY:
            self.logger.info (f"Not found next, first {local_scores[best_ind]}, second best {local_scores[second_best_ind]}")
            return "NONE"
        elif self.scaffold_graph.nodes[local_scores[best_ind][0]]['telomere'][0]:
            self.logger.info (f"Best {local_scores[best_ind]}, is good count but actually are in telomere!")            
            return "NONE"
        else:
            self.logger.info (f"Really best {local_scores[best_ind]}, second best {local_scores[second_best_ind]}")            
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
        self.logger.info ("Starting paths")
        self.logger.info (starting_paths)
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
                    self.logger.info (f"Extention {next_path_id} looks good but already used")
                    break
                if self.rc_path_id(self.findExtension(self.rc_path_id(next_path_id))) != cur_path_id:
                    self.logger.info (f"backward check failed for {next_path_id}")
                    break
                self.logger.info (f"Extending {cur_path_id} with {next_path_id}")

                #possibly not do so with paths of length one? They can be more successful in other direction
                cur_scaffold.append(next_path_id)
                nor_used_path_ids.add(next_path_id.strip('-+'))
                cur_path_id = next_path_id
            self.logger.info(f"scaffold {cur_scaffold}\n")
            res.append(cur_scaffold)
        total_scf = 0
        total_jumps = 0
        total_new_t2t = 0     
        final_paths = path_storage.PathStorage()   
        for nor_path_id in self.rukki_paths.getPathIds():
            if not (nor_path_id in nor_used_path_ids):
                final_paths.addPath(self.rukki_paths.getPathTsv(nor_path_id) , self.G)
        for scf in res:
            scf_path = []
            if len(scf) > 1:
                total_scf += 1
            total_jumps += len(scf) - 1
            largest_label = "NA"
            largest_id = scf[0].strip('-+')
            largest_len = 0
            cur_path_count = 0
            for or_path_id in scf:
                cur_path_count += 1
                nor_path_id = or_path_id.strip('-+')
                if self.rukki_paths.getLength(nor_path_id) > largest_len and self.rukki_paths.getLabel(nor_path_id) != "NA":
                    largest_len = self.rukki_paths.getLength(nor_path_id)
                    largest_id = nor_path_id
                    largest_label = self.rukki_paths.getLabel(nor_path_id)
                if or_path_id[-1] == "+":
                    scf_path.extend(self.rukki_paths.getPathById(nor_path_id))
                else:
                    scf_path.extend(sf.rc_path(self.rukki_paths.getPathById(nor_path_id)))
                if cur_path_count < len(scf):
                    scf_path.append("[N1000001N:scaffold]")
                    
            #debugging part
            for or_path_id in scf:
                nor_path_id = or_path_id.strip('-+')
                if nor_path_id in self.haploids and len(scf) > 1:
                    self.logger.info (f"scaffold part {nor_path_id} is haploid")
            for i in range (0, len(scf) - 1):
                if (scf[i].strip('-+'), scf[i+1].strip('-+')) in self.dangerous_swaps:
                    self.logger.warning (f"consecutive pair,  {scf[i]} {scf[i+1]} did signficant changes on hi-c counts based on {self.dangerous_swaps[(scf[i].strip('-+'), scf[i+1].strip('-+'))]}")
            path_str = "\t".join([largest_id, ",".join(scf_path), largest_label])
            final_paths.addPath(path_str, self.G)
            telo_end = False
            telo_start = False
            for tel_node in self.tel_nodes:
                if self.upd_G.has_edge(scf_path[-1], tel_node):
                    telo_end = True
                if self.upd_G.has_edge(tel_node, scf_path[0]):
                    telo_start = True
            if telo_end and telo_start:
                total_new_t2t += 1
            if len(scf) > 1:
                self.logger.info (f"SCAFFOLD {scf} {telo_start} {telo_end} ")
        self.logger.warning (f"Total scaffolds {total_scf} total jumps {total_jumps} new T2T {total_new_t2t}")
        self.outputScaffolds(final_paths)
        return res
    
    def outputScaffolds(self, final_paths):
        output_tsv = self.output_basename + ".tsv"
        output_gaf = self.output_basename + ".gaf"
        header = "name\tpath\tassignment"
        with open(output_tsv, "w") as tsv_file, open(output_gaf, "w") as gaf_file:
            tsv_file.write(header + "\n")
            gaf_file.write(header + "\n")
            for path_id in sorted(final_paths.getPathIds()):
                tsv_file.write(final_paths.getPathTsv(path_id) + "\n")
                gaf_file.write(final_paths.getPathGaf(path_id) + "\n")
        return
    #returns: dict {(start_id, end_id):[[start_pos1, end_pos1]]}. Coords not compressed!
    def get_connections(self, alignment_file):
        res = {}
        #A01660:39:HNYM7DSX3:1:1101:1696:29982   utig4-73        utig4-1056      1       16949880        78591191
        for line in open (alignment_file):
            arr = line.split()

            if not (arr[1], arr[2]) in res:
                res[(arr[1], arr[2])] = []
            next = [int(arr[4]), int(arr[5])]
            res[(arr[1], arr[2])].append(next)

            if not (arr[2], arr[1]) in res:
                res[(arr[2], arr[1])] = []
            res[(arr[2], arr[1])].append([int(arr[5]), int(arr[4])])  
        return res
    
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
        paths = [self.rukki_paths.getPathById(path_ids[0]), self.rukki_paths.getPathById(path_ids[1])]
        if total_score > 0:
            #whether we need to switch orientation for path number i
            for i in range (0, 2):
                #Do we need length check here at all?
                if self.rukki_paths.getLength(path_ids[i]) <= self.MAX_REORDERING_LENGTH:
                    correct_or = ""
                    if self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1]:
                        correct_or = correct_orientations[i]
                    elif self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0]:
                        correct_or = self.rc_orientation(correct_orientations[i])
                    if correct_or != "":            
                        pair_or = ["",""]
                        pair_or[i] = correct_or                        
                        self.logger.debug (f"tuning pair {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                        for second_or in ('-', '+'):
                            correct_pair = pair_or.copy()                       
                            correct_pair[1 - i] = second_or
                            incorrect_pair = correct_pair.copy()
                            incorrect_pair[i] = self.rc_orientation(correct_or)
                            correct_pair = "".join(correct_pair)
                            incorrect_pair = "".join(incorrect_pair)                            
                            #just logging
                            if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS and scores[incorrect_pair]  > scores[correct_pair]:
                                self.dangerous_swaps[(path_ids[0], path_ids[1])] = "telomeres"

                                self.logger.debug(f"Dangerous telomeric tuning pair {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                                self.logger.debug(f"Dangerous moved {scores[incorrect_pair]} from {incorrect_pair} to {correct_pair}")    
                            self.logger.debug (f"moving {incorrect_pair} to {correct_pair}")
                            scores[correct_pair] += scores[incorrect_pair]
                            scores[incorrect_pair] = 0                            
                        self.logger.debug (f"telomeric tuned pair {path_ids}, scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}") 
        
            for i in range (0, 2):
                #Do we need length check here? Should it be same as for telomeric one
                if self.rukki_paths.getLength(path_ids[i]) <= self.MAX_REORDERING_LENGTH:                  
                    for fixed_orientation in ('-', '+'):
                        shortest_paths = {'-':1000000000, '+':1000000000}
                        min_cutoff = min(self.CLOSE_IN_GRAPH, self.rukki_paths.getLength(path_ids[i]) / 4)                        
                        max_cutoff = self.rukki_paths.getLength(path_ids[i]) * 3 / 4
                        for orient in ('-', '+'):
                            to_check = paths.copy()
                            if fixed_orientation == "-":
                                to_check[1 - i] = sf.rc_path(paths[1-i])
                            if orient == '-':
                                to_check[i] = sf.rc_path(paths[i])
                            shortest_paths[orient] = self.pathDist(to_check[0], to_check[1], True)
                            self.logger.debug(f"Checking dists {to_check} index {i} dist {shortest_paths[orient]} cutoffs {min_cutoff} {max_cutoff}")

                        if shortest_paths['-'] < min_cutoff and shortest_paths['+'] > max_cutoff:
                            correct_or = "-"    
                        elif shortest_paths['+'] < min_cutoff and shortest_paths['-'] > max_cutoff:
                            correct_or = "+"
                        else:
                            correct_or = ""
                        if correct_or != "":
                            correct_pair = ["", ""]
                            incorrect_pair = ["", ""]
                            correct_pair[i] = correct_or
                            correct_pair[1 - i] = fixed_orientation
                            incorrect_pair[i] = self.rc_orientation(correct_or)
                            incorrect_pair[1 - i] = fixed_orientation   
                            correct_pair = "".join(correct_pair)
                            incorrect_pair = "".join(incorrect_pair) 
                            #just logging
                            if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS and scores[incorrect_pair] > scores[correct_pair]:
                                self.dangerous_swaps[(path_ids[0], path_ids[1])] = "connectivity"

                                self.logger.debug(f"Dangerous connectivity tuning pair {path_ids}, i {i} scores {scores}")                    
                                self.logger.debug(f"Dangerous moved {scores[incorrect_pair]} from {incorrect_pair} to {correct_pair}")    
                        
                            scores[correct_pair] += scores[incorrect_pair]
                            scores[incorrect_pair] = 0 
                            self.logger.debug (f"Connectivity tuned pair {path_ids}, scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}") 


        return scores

    #return scores for each of the orientations, ++, -+, +-, --,
    #orientation within the path, it is not changing!
    def getNodePairConnections(self, pair, orientations, connections, shift_before, shift_after, lens):    
        #This is scores for PATH orientation
        scores = {"++":0, "-+":0, "+-":0, "--":0, "middle":0}    
        filtered = 0
        not_filtered = 0
        if self.matchGraph.has_edge(pair[0], pair[1]):
            intervals = self.matchGraph.edges[pair[0], pair[1]]['intervals']
        else:
            intervals = [[],[]]
        for conn in connections[pair]:        
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
                    if (dists_to_end[i][j]< ScaffoldGraph.NEAR_PATH_END):
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
        self.logger.debug (f"Scores for {pair} {scores} {orientations} filtered/not_filtered {filtered} {not_filtered}")
        self.logger.debug (f"Shifts {shift_before} {shift_after}")
        self.logger.debug (f"Ignored intervals {intervals}")
        return scores


    #return scores for each of the orientations, ++, -+, +-, --,
    def getPathPairConnections(self, path_ids, connections, lens):
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
            if not (nor_f in lens) or lens[nor_f] < ScaffoldGraph.SHORT_INGORED_NODE:
                continue
            if not first in self.multiplicities or self.multiplicities[first] > 1:
                continue
            for second in paths[1]:
                nor_s = second.strip('-+')
                if not (nor_s in lens) or lens[nor_s] < ScaffoldGraph.SHORT_INGORED_NODE:
                    continue
                if not second in self.multiplicities or self.multiplicities[second] > 1:
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
                scores_pair = self.getNodePairConnections(pair, orientations, connections, shift_before, shift_after, lens)
                for key in scores_pair:
                    scores[key] += scores_pair[key] 
        self.logger.debug (f"Pathscores for {path_ids} {scores}")
        scores = self.fixOrientation(path_ids, scores)
        return scores
