#!/usr/bin/env python3

import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions as gf
import logging
from scaffolding import logger_wrap, match_graph, path_storage

class ReferencePosition:
    def __init__(self, name_q, name_r, ref_start, ref_end, query_len, orientation):
        self.name_q = name_q
        self.name_r = name_r
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.average_pos = (ref_end + ref_start)/2
        self.query_len = query_len
        self.orientation = orientation

#TODO: or inherit from nx.Digraph??
class ScaffoldGraph:
    LONG_HAPLOID_CUTOFF = 5000000
    #should not be actuallly used, since they'll be reordered with MAX_REORDERING_LENGTH check
    SHORT_TEL_CUTOFF = 5000000

    #Should it be ignored at all? 
    MAX_REORDERING_LENGTH = 30000000

    #If not near end, hi-c link is ignored
    NEAR_PATH_END = 5000000
    #<MIN_LINKS - ignore even clear majority
    MIN_LINKS = 10
    #Ratio from best to second best. Possibly should be increased.
    CLEAR_MAJORITY = 1.5   
    SHORT_INGORED_NODE = 20000
    #ignore shorter paths, TODO change for iterative runs
    MIN_PATH_TO_SCAFFOLD = 200000

    #Distance is defined with respect to homologous paths to allow "gap jumping"
    #if one of orientation is relatively close in graph(<min(1/4*path_length, CLOSE_IN_GRAPH) and other is really far (>3/4 of the length), we move all connections from far one to close
    CLOSE_IN_GRAPH = 500000
    #If paths are closer than CLOSE_IN_GRAPH, we significantly increase scores. Should be reconsidered when lots of gaps present
    #can be asymetric because of the 1/4 path_length rule, possibly should reconsider it
    CONNECTIVITY_MULTIPLICATIVE_BONUS = 2


    #default values for MatchGraph construction
    MATCHGRAPH_LONG_NODE = 500000
    MATCHGRAPH_MIN_ALIGNMENT = 100000
    

    #reference-based/haplotype-reference based params
    #shorter homologies will be ignored, possibly can be tuned in mashmap options?
    #TODO: shouldn't be always same as MIN_PATH_TO_SCAFFOLD?
    MIN_HOMOLOGY_REF = 100000
    #Quite conservative here, best homology contig should be twice longer than second best
    RATIO_HOMOLOGY_REF = 2.0
    #If best alignment is not covered by homologous intervals good enough - ignore 
    #TODO not sure whether we need it
    LEN_FRAC_REF = 0.6
    #Giving ref bonus for neighboring paths if dist is smaller than
    ABSOLUTE_ALLOWED_REFERENCE_JUMP = 5000000
    #consecutive by reference can not overlap longer than this
    ABSOLUTE_ALLOWED_REFERENCE_OVERLAP = MIN_PATH_TO_SCAFFOLD
    #If there is a path between that is longer than this * jump_in_reference_coords - then paths we are checking are considered to be not consistent
    INTERMEDIATE_PATH_FRACTION = 0.7
    
    #Consequtive paths scores are increased by this factor. 
    #TODO Possibly should be some high absolute constant too?
    REFERENCE_MULTIPLICATIVE_BONUS = 4
    


    def __init__(self, rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, uncompressed_fasta, path_mashmap, logger):
        self.logger = logger_wrap.UpdatedAdapter(logger, self.__class__.__name__)
        self.rukki_paths = rukki_paths

        self.multiplicities = rukki_paths.getEdgeMultiplicities()
        #upd_G - graph with telomeric nodes
        self.tel_nodes, self.upd_G = gf.get_telomeric_nodes(telomere_locations_file, G)
        self.logger.debug("Telomeric nodes")
        self.logger.debug(self.tel_nodes)
        self.output_basename = "scaff_rukki.paths"

        #Used for scaffolding starts and debug, 
        self.to_scaff = self.get_paths_to_scaff(ScaffoldGraph.MIN_PATH_TO_SCAFFOLD)
        self.hic_alignment_file = hic_alignment_file

        #TODO: duplicated... self.matchGraph should not be refered directly
        self.mg = match_graph.MatchGraph(matches_file, G, -239239239, ScaffoldGraph.MATCHGRAPH_LONG_NODE, ScaffoldGraph.MATCHGRAPH_MIN_ALIGNMENT, logger)
        self.matchGraph = self.mg.getMatchGraph()
        
        #debug purposes
        self.dangerous_swaps = {}

        self.uncompressed_lens = gf.get_lengths(uncompressed_fasta)
        self.compressed_lens = {}
        for node in self.upd_G.nodes:
            self.compressed_lens[node] = self.upd_G.nodes[node]['length']      
            self.compressed_lens[node.strip("-+")] = self.upd_G.nodes[node]['length']  

        #positions on refence or other haplotype  
        self.reference_positions = {}
        self.assigned_reference = {}
        self.getPathPositions(path_mashmap)
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

#TODO: move all rc_<smth> somewhere, not right place 
    def rc_orientation(self, c):
        if c == "+":
            return "-"
        if c == "-":
            return "+"
        return c    

    def rc_path_id(self, path_id):
        return path_id[:-1] + self.rc_orientation(path_id[-1])
    
    def rc_path(self, path):
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
                #DEBUG ONLY
                if path_len > 2000000:
                    self.logger.info(f"Found haploid path {nor_path_id} with homology {total_hom} and len {path_len} ")
        return haploids

        
    def getClosestTelomere(self, path, direction):
        #From telomere to rc(path_end)
        if direction == '+':
            path = self.rc_path(path)
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
            #self.logger.debug (f"Checking nodepair dist from {node} to {path_node}")
            if path_node in self.upd_G.nodes:
                if path_node == node:
                    closest = 0
                else:
                    if path_node in self.dists[node]:
                        closest = min(closest, self.dists[node][path_node] + add_dist - self.compressed_lens[node.strip("-+")] - self.compressed_lens[path_node.strip("-+")])
                    if check_homologous:
                        #self.logger.debug(f"Checking homologous to node {path_node}: {self.homologousOrNodes(path_node)}")
                        for hom_node in self.homologousOrNodes(path_node):
                            #self.logger.debug (f"Checking homologous nodepair dist from {node} to {hom_node}")
                            if hom_node == node:
                                closest = 0
                            elif hom_node in self.dists[node]:
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
                        #self.logger.debug (f"Checking homologous dist from {hom_node} to {path_to} add_dist {add_dist}")
                        closest = min(closest, (self.nodeToPathDist(hom_node, path_to, check_homologous) + add_dist))
                add_dist += self.compressed_lens[node.strip("-+")]
        return closest

    def getPathPositions(self, path_mashmap):
        hom_storage = match_graph.HomologyStorage(self.logger, path_mashmap, ScaffoldGraph.MIN_HOMOLOGY_REF)
        for path_id in hom_storage.homologies:
            all_refs = []
            if self.rukki_paths.getLength(path_id) < ScaffoldGraph.MIN_PATH_TO_SCAFFOLD:
                continue
            for ref_id in hom_storage.homologies[path_id]:
                all_refs.append([hom_storage.homologies[path_id][ref_id].getCoveredLen(), ref_id])
            all_refs.sort(reverse=True)
            if len(all_refs) == 1 or all_refs[0][0] > ScaffoldGraph.RATIO_HOMOLOGY_REF * all_refs[1][0]:
                best_ref = all_refs[0][1]
                self.logger.debug(f"Clearly best homology for {path_id} is {best_ref} with size{all_refs[0][0]}")
                if all_refs[0][0] > ScaffoldGraph.LEN_FRAC_REF * hom_storage.getLength(path_id):
                    self.logger.debug(f"Best homology for {path_id} is {best_ref} with size{all_refs[0][0]}")
                    if not (best_ref+'+') in self.reference_positions:
                        self.reference_positions[best_ref+ "+"] = []
                        self.reference_positions[best_ref+ "-"] = []
                    hom_info = hom_storage.homologies[path_id][best_ref] 
                    self.assigned_reference[path_id + '+'] = best_ref + hom_info.orientation
                    self.assigned_reference[path_id + '-'] = best_ref + self.rc_orientation(hom_info.orientation)
#                    self.assigned_reference[path_id] = best_ref
#                    self.reference_positions[best_ref].append(ReferencePosition(path_id, best_ref, hom_info.largest_interval_center[1], hom_storage.getLength(path_id), hom_info.orientation))
                    
                    self.reference_positions[best_ref + "+"].append(ReferencePosition(path_id + hom_info.orientation, best_ref + "+", hom_info.approximate_positions[1][0],hom_info.approximate_positions[1][1], hom_storage.getLength(path_id), hom_info.orientation))
                    self.reference_positions[best_ref + "-"].append(ReferencePosition(path_id + self.rc_orientation(hom_info.orientation), best_ref + '-', hom_storage.getLength(best_ref) - hom_info.approximate_positions[1][1], hom_storage.getLength(best_ref) - hom_info.approximate_positions[1][0], hom_storage.getLength(path_id), self.rc_orientation(hom_info.orientation)))

                else:
                    self.logger.debug(f"Best homology for {path_id} is {best_ref} not covered enough frac {all_refs[0][0] / hom_storage.getLength(path_id)}")
        for ref in self.reference_positions:
            self.logger.debug(f"Reference positions for {ref}")
            self.reference_positions[ref] = sorted(self.reference_positions[ref], key=lambda x: x.average_pos)
            for pos in self.reference_positions[ref]:
                self.logger.debug(f"{pos.name_q} {pos.name_r} {pos.average_pos} {pos.ref_start} {pos.ref_end} {pos.query_len} {pos.orientation}")          


    def isNextByRef(self, from_path_id, to_path_id):
        if (not (from_path_id in self.assigned_reference) or not (to_path_id in self.assigned_reference)):
            return False
        if self.assigned_reference[from_path_id] != self.assigned_reference[to_path_id]:
            return False
        aligned = self.reference_positions[self.assigned_reference[from_path_id]]    
        len_aligned = len(aligned)
        for i in range(0, len_aligned - 1):
            if aligned[i].name_q == from_path_id:
                #can be not consequtive if reference is used - both haplos are aligned
                for j in range (i +1, len_aligned):
                    if aligned[j].name_q == to_path_id:
                        prev_end = aligned[i].ref_end
                        next_start = aligned[j].ref_start
                        inconsistency_len = (next_start - prev_end)
                        if inconsistency_len > ScaffoldGraph.ABSOLUTE_ALLOWED_REFERENCE_JUMP:
                            return False
                        #TODO: possibly should be tuned, since positions are not exact. Or more accurate position computation....
                        if inconsistency_len < -1 * ScaffoldGraph.ABSOLUTE_ALLOWED_REFERENCE_OVERLAP:
                            return False
                        if inconsistency_len > 0:
                            for mid in range (i+1, j):
                                #there is some path in between of i and j, they are not consistent
                                if aligned[mid].query_len * ScaffoldGraph.INTERMEDIATE_PATH_FRACTION < inconsistency_len:
                                    return False
                        return True
                        
        return False

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
        final_paths = path_storage.PathStorage(self.upd_G)   
        for nor_path_id in self.rukki_paths.getPathIds():
            if not (nor_path_id in nor_used_path_ids):
                final_paths.addPath(self.rukki_paths.getPathTsv(nor_path_id))
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
                    scf_path.extend(self.rc_path(self.rukki_paths.getPathById(nor_path_id)))
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
            final_paths.addPath(path_str)
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
                                to_check[1 - i] = self.rc_path(paths[1-i])
                            if orient == '-':
                                to_check[i] = self.rc_path(paths[i])
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
                            scores[correct_pair] *= self.CONNECTIVITY_MULTIPLICATIVE_BONUS
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
        #self.logger.debug (f"Scores for {pair} {scores} {orientations} filtered/not_filtered {filtered} {not_filtered}")
        #self.logger.debug (f"Shifts {shift_before} {shift_after}")
        #self.logger.debug (f"Ignored intervals {intervals}")
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
                #not counting homozygous nodes, they can be really large
                if or_node.strip('-+') in lens and or_node in self.multiplicities and self.multiplicities[or_node] == 1:
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
        for orientation in scores:
            if len(orientation) == 2:
                if self.isNextByRef(path_ids[0] + orientation[0], path_ids[1] + orientation[1]):
                    self.logger.debug (f"Reference connection found! {path_ids} {orientation}")
                    scores[orientation] *= self.REFERENCE_MULTIPLICATIVE_BONUS
        return scores
    
    #returns dict, {id:[present_start_relo, present_end_telo]}
    def get_paths_to_scaff(self, long_enough):
        res = {}
        for id in self.rukki_paths.paths:
            total_l = self.rukki_paths.path_lengths[id]
            tel_start = False
            tel_end = False    
            for telomere in self.tel_nodes:
                if self.upd_G.has_edge(telomere, self.rukki_paths.paths[id][0]):
                    tel_start = True
                if self.upd_G.has_edge(self.rukki_paths.paths[id][-1], telomere):
                    tel_end = True
            if tel_start and tel_end:
                continue
            #all long enough AND containing telomere
            if total_l > long_enough:
                res[id] = [tel_start, tel_end]
                #print (f"will use path {paths.paths[id]} {tel_start} {tel_end}")
        return res