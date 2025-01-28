#!/usr/bin/env python3

from collections import namedtuple
import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions as gf
import logging

#from pympler import asizeof

#from memory_profiler import profile
import pysam
#TODO: remove clother to release? 
import psutil
from scaffolding import logger_wrap, match_graph, path_storage
#from line_profiler import LineProfiler

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
    #to avoid storing doubles as multimpaper weight will multiply by this and round, possibly to remove
    INT_NORMALIZATION = 100


    #TODO possibly change the definition - whether we do not have enough similarity at all, or do not have one2one similarity 
    #Want this to be longer than longest distal bit, heterogametic sex chrs only
    LONG_HAPLOID_CUTOFF = 6000000
    #should not be actuallly used, since they'll be reordered with MAX_REORDERING_LENGTH check
    #not in use except of assertion
    SHORT_TEL_CUTOFF = 5000000

    #Should it be ignored at all? 
    MAX_REORDERING_LENGTH = 30000000

    #If not near end, hi-c link is ignored
    NEAR_PATH_END = 5000000
    #<MIN_LINKS - ignore even clear majority
    MIN_LINKS = 10 
    #Ratio from best to second best. Possibly should be increased.
    CLEAR_MAJORITY = 1.5   
    #TODO: sinice multiplicity > 1 is not used, do we really need to ignore short ones?
    SHORT_INGORED_NODE = 20000
    #ignore shorter paths, TODO change for iterative runs
    MIN_PATH_TO_SCAFFOLD = 200000

    #Consts connecting two remaining telo-containing paths in component
    MIN_EXPECTED_T2T = 1000000

    DEFAULT_GAP_SIZE = 100000
    #Distance is defined with respect to homologous paths to allow "gap jumping"
    #if one of orientation is relatively close in graph(<min(1/4*path_length, CLOSE_IN_GRAPH) and other is really far (>3/4 of the length), we move all connections from far one to close
    CLOSE_IN_GRAPH = 500000
    #If paths are closer than CLOSE_IN_GRAPH, we significantly increase scores. Should be reconsidered when lots of gaps present

    #can be asymetric because of the 1/4 path_length rule, possibly should reconsider it
   #TODO: reconsidered to be smaller than min of two paths, check whether anything go wrong

   


    #default values for MatchGraph construction
    MATCHGRAPH_LONG_NODE = 500000
    MATCHGRAPH_MIN_ALIGNMENT = 100000
    

    #reference-based/haplotype-reference based params
    #shorter homologies will be ignored, possibly can be tuned in mashmap options?
    #TODO: shouldn't be always same as MATCHGRAPH_MIN_ALIGNMENT? Rename REFERENCE_MIN_ALIGNMENT? Likely decrease to 50K
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
    #Without it paths from alternative haplotype can interfere 
    INTERMEDIATE_PATH_FRACTION = 0.7
    
    #Consequtive paths scores are increased by this factor. 
    #TODO Possibly should be some high absolute constant too?
    REFERENCE_MULTIPLICATIVE_BONUS = 5
    #efficiently was 4 in lots of cases because applied twice
    #TODO recheck
    CONNECTIVITY_MULTIPLICATIVE_BONUS = 3

    #Just too long/far
    TOO_FAR = 1000000000

    
    #efficiantly it is so because of bwa behaviour on XA tags but not used directly
    MAX_ALIGNMENTS = 6

    APPROXIMATE_COORDS = 40000
    APPROXIMATE_COORDS_HALF = APPROXIMATE_COORDS/2

    #to check whether paths have similar lengths, in cheating t2t connection
    SIMILAR_LEN_FRAC = 0.7

    #Constants to switch off distance counting for huge graphs
    MAX_GRAPH_FOR_DISTANCES = 1000000
    MAX_COMPONENT_FOR_DISTANCES = 50000

    def __init__(self, rukki_paths, telomere_locations_file, hic_alignment_file, matches_file, G, uncompressed_fasta, path_mashmap, porec, logger):
        self.logger = logger_wrap.UpdatedAdapter(logger, self.__class__.__name__)
        self.rukki_paths = rukki_paths        
        self.uncompressed_lens = gf.get_lengths(uncompressed_fasta)
        self.multiplicities = rukki_paths.getEdgeMultiplicities()
        interesting_nodes = self.getInterestingNodes(self.uncompressed_lens)
        self.logger.info(f"Total nodes {len(G.nodes)} interesting nodes {len(interesting_nodes)}")
        
        #upd_G - graph with telomeric nodes
        self.tel_nodes, self.upd_G = gf.get_telomeric_nodes(telomere_locations_file, G)
        self.logger.debug("Telomeric nodes")
        self.logger.debug(self.tel_nodes)

        self.match_graph = match_graph.MatchGraph(matches_file, G, -239239239, ScaffoldGraph.MATCHGRAPH_LONG_NODE, ScaffoldGraph.MATCHGRAPH_MIN_ALIGNMENT, logger)
        
        if porec == True:
            self.logger.info ("Loading pore-c alignments")
            all_connections, unique_connections = self.get_connections_porec(hic_alignment_file, True)
        else:
            self.logger.info ("Loading Hi-C alignments")
            all_connections, unique_connections = self.get_connections_bam(hic_alignment_file, True) 

        
        self.all_connections = all_connections
        self.unique_connections = unique_connections
        self.all_connections_nodemap = self.getConnectionsNodemap(all_connections)
        self.unique_connections_nodemap = self.getConnectionsNodemap(unique_connections)

        self.output_basename = "scaff_rukki.paths"
        telomeric_ends = self.getTelomericEnds()
        self.hic_alignment_file = hic_alignment_file

        
        #debug purposes
        self.dangerous_swaps = {}

        self.compressed_lens = {}
        for node in self.upd_G.nodes:
            self.compressed_lens[node] = self.upd_G.nodes[node]['length']      
            self.compressed_lens[node.strip("-+")] = self.upd_G.nodes[node]['length']  

        #positions on refence or other haplotype  
        self.reference_positions = {}
        self.assigned_reference = {}
        self.getPathPositions(path_mashmap)
        self.G = G
  
        self.scaffold_graph = nx.DiGraph()

        for id in rukki_paths.getPathIds():
            for dir in ('-', '+'):
                or_id = id + dir
                if dir == '+':
                    tels_ends = telomeric_ends[id]
                else:
                    tels_ends = [telomeric_ends[id][1], telomeric_ends[id][0]]
                self.scaffold_graph.add_node(or_id, telomere = tels_ends)    
        #possilby unefficient but whelp
        scores = {}
        unique_scores = {}
        self.dists = {}
        if G.number_of_nodes() > ScaffoldGraph.MAX_GRAPH_FOR_DISTANCES:
            self.logger.info(f"Graph is too big {G.number_of_nodes()}, not calculating pairwise distances")
        else:
            max_comp_size = len(max(nx.weakly_connected_components(G), key=len))    
            if max_comp_size > ScaffoldGraph.MAX_COMPONENT_FOR_DISTANCES:
                self.logger.info(f"Biggest component is too big {max_comp_size}, not calculating pairwise distances")
            else:
                self.logger.info("Calculating pairwise distances for assembly graph nodes")
                self.dists = dict(nx.all_pairs_dijkstra_path_length(self.upd_G, weight=lambda u, v, d: self.upd_G.edges[u, v]['mid_length']))
                self.logger.info("Pairwise distances in assembly graph calculated")

        self.haploids = self.getHaploidPaths(self.rukki_paths)
        self.scores = {}
        self.unique_scores = {}

    
    def getConnectionsNodemap(self, connections):
        nodemap = {}
        for node_pair in connections.keys():
            if not node_pair[0] in nodemap:
                nodemap[node_pair[0]] = set()
            if not node_pair[1] in nodemap:
                nodemap[node_pair[1]] = set()          
            nodemap[node_pair[0]].add(node_pair[1])
            nodemap[node_pair[1]].add(node_pair[0])
        return nodemap
    
    def isHaploidPath(self, paths, nor_path_id):
        total_hom = 0
        for or_node in paths.getPathById(nor_path_id):
            nor_node = or_node.strip("-+")

            for next in self.match_graph.getHomologousNodes(nor_node, True):
                total_hom += self.match_graph.getEdgeAttribute(nor_node, next, 'homology_len')
        path_len = paths.getLength(nor_path_id)
        #TODO: should exclude homozygous nodes here
        if total_hom * 2 < path_len:            
            #DEBUG ONLY                
            if path_len > 2000000:
                self.logger.info(f"Found haploid path {nor_path_id} with homology {total_hom} and len {path_len} ")
            return True
        else:
            return False

    def getHaploidPaths(self, paths):
        haploids = set()
        nodes_to_path_ids = {}        
        for nor_path_id in paths.getPathIds():
            for or_node in paths.getPathById(nor_path_id):
                nor_node = or_node.strip("-+")
                if not (nor_node in nodes_to_path_ids):
                    nodes_to_path_ids[nor_node] = []
                nodes_to_path_ids[nor_node].append(nor_path_id)
        #possibly we'll need that graph but not now
        #homGraph = nx.Graph()
        for nor_path_id in paths.getPathIds():
            if self.isHaploidPath(paths, nor_path_id):
                haploids.add(nor_path_id)
        return haploids
    
    #TODO a lot of graph-only code should be moved in or_graph wrapper class
    #we store double dists from middle to middle, transforming to dist from end to    
    def orNodeDist(self, from_node, to_node):
        res = ScaffoldGraph.TOO_FAR
        if from_node in self.dists:
            if to_node in self.dists[from_node]:
                res = (self.dists[from_node][to_node] - self.compressed_lens[gf.nor_node(from_node)] - self.compressed_lens[gf.nor_node(to_node)])/2
        return res
    
    #Dist from node end to path. Allowed to go not in the first path node, but then additional length in path is added
    #Optionally allowing to use homologous nodes (to improve in gaps)
    def nodeToPathDist(self, node, path, check_homologous):
        closest = ScaffoldGraph.TOO_FAR 
        add_dist = 0
        for path_node in path:
            #self.logger.debug (f"Checking nodepair dist from {node} to {path_node}")
            if path_node in self.upd_G.nodes:
                #TODO: likely this is not needed
#                if path_node == node:
#                    closest = 0
#                else:                
                closest = min(closest, add_dist + self.orNodeDist(node, path_node))
                if check_homologous:
                    #self.logger.debug(f"Checking homologous to node {path_node}: {self.homologousOrNodes(path_node)}")
                    for hom_node in self.match_graph.getHomologousOrNodes(path_node, True):
                        #self.logger.debug (f"Checking homologous nodepair dist from {node} to {hom_node}")
                        if hom_node == node:
                            closest = 0                        
                        closest = min(closest, add_dist + self.orNodeDist(node, hom_node))
                add_dist += self.compressed_lens[gf.nor_node(path_node)]
        return closest

    #not to be used directly, with IDs we presave dists in pathsStorage
    def pathDist(self, path_from:list, path_to:list, check_homologous:bool):
        closest = ScaffoldGraph.TOO_FAR
        add_dist = 0
        for node in reversed(path_from):
            if node in self.upd_G.nodes:
                closest = min(closest, self.nodeToPathDist(node, path_to, check_homologous) + add_dist)
                if check_homologous:
                    for hom_node in self.match_graph.getHomologousOrNodes(node, True) :
                        #self.logger.debug (f"Checking homologous dist from {hom_node} to {path_to} add_dist {add_dist}")
                        closest = min(closest, (self.nodeToPathDist(hom_node, path_to, check_homologous) + add_dist))
                add_dist += self.compressed_lens[node.strip("-+")]
        return closest
    
    #Optionally allowing to use homologous nodes (to improve in gaps)
    def orPathIdDist(self, or_path_id_from:str, or_path_id_to:str, paths:path_storage.PathStorage, check_homologous:bool):
        precounted = paths.getStoredDist(or_path_id_from, or_path_id_to, check_homologous)
        if precounted != -1:
            return precounted
        ids = [or_path_id_from, or_path_id_to]
        to_dists = []
        for k in range (0, 2):
            if ids[k][-1] == "-":
                to_dists.append(gf.rc_path(paths.getPathById(ids[k][:-1])))
            else:
                to_dists.append(paths.getPathById(ids[k][:-1]))
        dist = self.pathDist(to_dists[0], to_dists[1], check_homologous)
        paths.storeDist(or_path_id_from, or_path_id_to, check_homologous, dist)
        return dist

    def norPathIdDist(self, nor_path_id_from:str, nor_path_id_to:str, paths:path_storage.PathStorage, check_homologous:bool):
        ids = [nor_path_id_from, nor_path_id_to]        
        res = ScaffoldGraph.TOO_FAR
        for orients in ('--', '++', '-+', '+-'):
            to_dists = []
            for k in range (0, 2):                
                to_dists.append(ids[k] + orients[k])
            res = min(res, self.orPathIdDist(to_dists[0], to_dists[1], paths, check_homologous))
        return res


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
                    self.assigned_reference[path_id + '-'] = best_ref + gf.rc_orientation(hom_info.orientation)
#                    self.assigned_reference[path_id] = best_ref
#                    self.reference_positions[best_ref].append(ReferencePosition(path_id, best_ref, hom_info.largest_interval_center[1], hom_storage.getLength(path_id), hom_info.orientation))
                    #TODO: possibly unaligned sequences after/before approximate interavals should be counted?
                    self.reference_positions[best_ref + "+"].append(ReferencePosition(path_id + hom_info.orientation, best_ref + "+", hom_info.approximate_positions[1][0],hom_info.approximate_positions[1][1], hom_storage.getLength(path_id), hom_info.orientation))
                    self.reference_positions[best_ref + "-"].append(ReferencePosition(path_id + gf.rc_orientation(hom_info.orientation), best_ref + '-', hom_storage.getLength(best_ref) - hom_info.approximate_positions[1][1], hom_storage.getLength(best_ref) - hom_info.approximate_positions[1][0], hom_storage.getLength(path_id), gf.rc_orientation(hom_info.orientation)))

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
            if self.orPathIdDist(from_path_id, to_path_id, self.rukki_paths, True) > ScaffoldGraph.NEAR_PATH_END:
                self.logger.warning(f"Banning distant link from long haploid {from_path_id} to diploid {to_path_id}")
                return True
            else:
                self.logger.warning(f"Allowing link from long haploid {from_path_id} to diploid {to_path_id}, close in graph")
        if nor_to_path_id in self.haploids and self.rukki_paths.getLength(nor_to_path_id) > ScaffoldGraph.LONG_HAPLOID_CUTOFF and not (nor_from_path_id in self.haploids):
            if self.orPathIdDist(from_path_id, to_path_id, self.rukki_paths, True) > ScaffoldGraph.NEAR_PATH_END:
                self.logger.warning(f"Banning distant link from diploid {from_path_id} to long haploid {to_path_id}")
                return True
            else:
                self.logger.warning(f"Allowing link from diploid {from_path_id} to long haploid {to_path_id}, close in graph")
        #relatively short fragments with telomere are special case, we may fail to detect orientation there but belive in telomere.
        if self.rukki_paths.getLength(nor_to_path_id) <= ScaffoldGraph.SHORT_TEL_CUTOFF and self.scaffold_graph.nodes[to_path_id]['telomere'][0]:
            self.logger.error(f"Banning link from {from_path_id} into short telomere, should not happen {to_path_id}")
            return True
        return False
    

    def getAnyLinkPaths(self, nor_path_id, path_storage, connection_map):
        connected_nodes = set()
        for node in path_storage.getPathById(nor_path_id):
            nor_node = gf.nor_node(node)
            if nor_node in connection_map:
                connected_nodes.update(connection_map[nor_node])
        connected_paths = set()
        for nor_node in connected_nodes:
            connected_paths.update(path_storage.getPathsFromNode(nor_node))
        return connected_paths


    def getScores(self, or_path_id, path_storage, type):
        res = []
        if type == "weight":
            cur_scores = self.scores
            connections = self.all_connections
            connection_nodemap = self.all_connections_nodemap
        elif type == "unique_weight":
            cur_scores = self.unique_scores
            connections = self.unique_connections
            connection_nodemap = self.unique_connections_nodemap

        else:
            self.logger.error(f"Unknown type {type}")
            return res
        nor_path_id = gf.nor_path_id(or_path_id)
        precounted = True

        homologous_paths = self.match_graph.getHomologousPaths(path_storage, nor_path_id)
        #Counting only necessary scores, cashing in global dict 
        if not (or_path_id in cur_scores):
            cur_scores[or_path_id] = {}
            cur_scores[gf.rc_path_id(or_path_id)] = {}
            precounted = False
        connected_nor_paths = self.getAnyLinkPaths(nor_path_id, path_storage, connection_nodemap)

        for next_nor_path in connected_nor_paths:
            if next_nor_path in homologous_paths or next_nor_path == nor_path_id:
                continue
            else:
                if not precounted:
                    nor_scores = self.getPathPairConnections([nor_path_id, next_nor_path], connections, self.uncompressed_lens)                
                    for from_dir in ('-', '+'):
                        for to_dir in ('-', '+'):
                            or_from_path_id = nor_path_id + from_dir
                            or_to_path_id = next_nor_path + to_dir
                            if nor_scores[from_dir + to_dir] > 0:
                                cur_scores[or_from_path_id][or_to_path_id] = nor_scores[from_dir + to_dir]
        
        for next_or_path in cur_scores[or_path_id].keys():
            res.append([next_or_path, cur_scores[or_path_id][next_or_path]])

        res.sort(key=lambda x: x[1], reverse=True)
        return res
    
    #Main logic is here!     
    # returns NONE/DONE/next_path_id   
    #type: weight/unique_weight
    def findExtension(self, cur_path_id, type):
        local_scores = []
        self.logger.info(f"Checking {cur_path_id}")
        if not (cur_path_id in self.scaffold_graph.nodes):
            return "NONE",0
        if self.scaffold_graph.nodes[cur_path_id]['telomere'][1]:
            self.logger.info (f"Stopped at the telomere")
            return "DONE",0
        local_scores = self.getScores(cur_path_id, self.rukki_paths, type)

        best_ind = -1
        second_best_ind = -1        
        for i in range (0, len(local_scores)):
            if (type != "unique_weight") or not self.forbiddenPair(cur_path_id, local_scores[i][0]):
                best_ind = i
                break
        for j in range (best_ind+1, len(local_scores)):
            if (type != "unique_weight") or not self.forbiddenPair(cur_path_id, local_scores[j][0]):
                second_best_ind = j
                break
        if len(local_scores) == 0:            
            self.logger.info (f"Nothing next")  
            return "NONE", 0
        elif local_scores[best_ind][1] <= ScaffoldGraph.MIN_LINKS:
            self.logger.info (f"very few links, best valid candidate {local_scores[best_ind]}")                 
            return "NONE", 0
        #not valid but waay best solution exists
        elif self.scaffold_graph.nodes[local_scores[best_ind][0]]['telomere'][0]:
            self.logger.info (f"Best {local_scores[best_ind]}, is good count but actually are in telomere!")            
            return "NONE", 0
        elif len(local_scores) == 1 or (len(local_scores) > 1 and  local_scores[second_best_ind][1] == 0):
            self.logger.info (f"Only one next, {local_scores[best_ind]}") 
            return local_scores[best_ind][0], ScaffoldGraph.CLEAR_MAJORITY*2                      
        else:
            self.logger.info (f"Really best {local_scores[best_ind]}, second best {local_scores[second_best_ind]}")            
            return local_scores[best_ind][0],local_scores[best_ind][1]/local_scores[second_best_ind][1]

    #made function to allow to use uniques
    #types: weight/unique_weight
    def findNextPath(self, current_path_id, nor_used_path_ids, type):
        next_path_id, ratio_fwd = self.findExtension(current_path_id, type)
        if next_path_id.strip('-+') in nor_used_path_ids:
            self.logger.info (f"Extention {next_path_id} looks good but already used")
            return "NONE"
        prev_path_id, ratio_bwd = self.findExtension(gf.rc_path_id(next_path_id), type)
        if gf.rc_path_id(prev_path_id) != current_path_id:
            self.logger.info (f"backward check failed for {next_path_id}")
            return "NONE"
        if ratio_fwd > ScaffoldGraph.CLEAR_MAJORITY and ratio_bwd > ScaffoldGraph.CLEAR_MAJORITY:
            self.logger.info (f"Extention {next_path_id} looks good")
            return next_path_id
        elif ratio_fwd * ratio_bwd > ScaffoldGraph.CLEAR_MAJORITY * ScaffoldGraph.CLEAR_MAJORITY * ScaffoldGraph.CLEAR_MAJORITY and type == "unique_weight":
            self.logger.info (f" Risky extension, not clear majority for one direction but very clear for other {current_path_id} {next_path_id} {ratio_fwd} {ratio_bwd}")
            return next_path_id
        return "NONE" 
        
    def generateScaffolds(self):
        res = []
        #will grow to right these paths in length order 
        tel_starting_paths = []        
        middle_paths = []
        #to avoid outputing same path twice
        nor_used_path_ids = set()
        for from_path_id in self.rukki_paths.getPathIds():
            if self.rukki_paths.getLength(from_path_id) < ScaffoldGraph.MIN_PATH_TO_SCAFFOLD:
                continue
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
        #from ID to amount of scaffolds
        scaff_sizes = {}
        for or_from_path_id in starting_paths:
            if or_from_path_id.strip('-+') in nor_used_path_ids:
                continue            
            cur_scaffold = [or_from_path_id]
            cur_path_id = or_from_path_id            
            nor_used_path_ids.add(or_from_path_id.strip('-+'))
            #lets add bidirectional expansion
            tried_reverse = False
            while True:
                next_path_id = self.findNextPath(cur_path_id, nor_used_path_ids, "weight")
                if next_path_id == "NONE":
                    self.logger.info (f"Failed to find regular extension for {cur_path_id}, trying unique")
                    next_path_id = self.findNextPath(cur_path_id, nor_used_path_ids, "unique_weight")

                if next_path_id == "DONE":
                    self.logger.info ("All done, stopped at telomere")
                    if tried_reverse:
                        break
                    else:
                        self.logger.info (f"Reversing {cur_scaffold}")
                        cur_scaffold = [gf.rc_path_id(x) for x in cur_scaffold]
                        cur_scaffold.reverse()
                        self.logger.info (f"Reversed {cur_scaffold}")
                        cur_path_id = cur_scaffold[-1]
                        tried_reverse = True
                        continue
                elif next_path_id == "NONE":
                    self.logger.info ("Failed to find extension, stopping")
                    if tried_reverse:
                        break
                    else:
                        self.logger.info (f"Reversing {cur_scaffold}")
                        cur_scaffold = [gf.rc_path_id(x) for x in cur_scaffold]
                        cur_scaffold.reverse()
                        self.logger.info (f"Reversed {cur_scaffold}")
                        cur_path_id = cur_scaffold[-1]
                        tried_reverse = True
                        continue
                    
                else:
                    hom_before = False
                    nor_new_path_id = gf.nor_path_id(next_path_id)
                    homologous_paths = self.match_graph.getHomologousPaths(self.rukki_paths, gf.nor_path_id(nor_new_path_id))
                    for prev_path_id in cur_scaffold:
                        nor_prev_path_id = gf.nor_path_id(prev_path_id)
                        if nor_prev_path_id in homologous_paths:
#                        self.match_graph.isHomologousPath([self.rukki_paths.getPathById(nor_new_path_id), self.rukki_paths.getPathById(nor_prev_path_id)], [self.rukki_paths.getLength(nor_new_path_id), self.rukki_paths.getLength(nor_prev_path_id)]):
                            #TODO: really not normal if we see that best extension is homologous to some path in scaffold, deserves investigation
                            self.logger.warning(f"Trying to extend, had homologous path in same scaffold before! {nor_new_path_id} {nor_prev_path_id}")
                            hom_before = True
                    if hom_before:
                        self.logger.info (f"Homologous paths before in scaffold, not extending {cur_path_id} with {next_path_id}")
                        if tried_reverse:
                            break
                        else:
                            cur_scaffold = [gf.rc_path_id(x) for x in cur_scaffold]
                            cur_scaffold.reverse()
                            cur_path_id = cur_scaffold[-1]
                            tried_reverse = True
                            continue                                            
                    self.logger.info (f"Extending {cur_path_id} with {next_path_id}")

                    #possibly not do so with paths of length one? They can be more successful in other direction
                    cur_scaffold.append(next_path_id)
                    nor_used_path_ids.add(next_path_id.strip('-+'))
                    cur_path_id = next_path_id
            #Reversing for comparison with previous runs only
            cur_scaffold = [gf.rc_path_id(x) for x in cur_scaffold]
            cur_scaffold.reverse()
            self.logger.info(f"scaffold {cur_scaffold}\n")
            if len(cur_scaffold) > 1:
                res.append(cur_scaffold)
            else:
                # not scaffolded, do not want it here because of rephasing
                nor_used_path_ids.remove(gf.nor_node(cur_path_id))

        total_scf = 0
        total_jumps = 0
        total_new_t2t = 0     
        final_paths = path_storage.PathStorage(self.upd_G)   
        scaffolded_rephased = []
        for scf in res:
            scaff_result = self.scaffoldFromOrPathIds(scf, self.rukki_paths)
            total_jumps += scaff_result[2]
            if scaff_result[2] > 0:
                total_scf += 1
            if scaff_result[1]:
                total_new_t2t += 1
            scaffolded_rephased.extend(scaff_result[3])
            final_paths.addPath(scaff_result[0])
        self.logger.info (f"Rephased scaffolds to other haplo, need to check homologous {scaffolded_rephased}")

        #this should be better done after completeConnectedComponent but it is very weird case if it matters         
        rephased_fix = self.getAdditionalRephase(scaffolded_rephased, nor_used_path_ids, self.rukki_paths)
        self.logger.info(f"Rephased fix {rephased_fix}")
        self.logger.info(f"used_ids {nor_used_path_ids}")
        
        for nor_path_id in self.rukki_paths.getPathIds():
            self.logger.debug(f"Checking {nor_path_id}")
            if not (nor_path_id in nor_used_path_ids):
                
                if nor_path_id in rephased_fix:
                    self.logger.debug(f"in reph fix {nor_path_id}")
                    final_paths.addPath(rephased_fix[nor_path_id])
                else:
                    final_paths.addPath(self.rukki_paths.getPathTsv(nor_path_id))

        self.logger.info("Starting last telomere tuning")
        component_tuned_paths = self.completeConnectedComponent(final_paths)
        telomere_cheating = len(final_paths.getPathIds()) - len(component_tuned_paths.getPathIds())
        self.fixHaploidNames(component_tuned_paths)        
        self.logger.warning (f"Total normal scaffolds {total_scf} last telomere tuned {telomere_cheating} total jumps {total_jumps + telomere_cheating} new T2T {total_new_t2t + telomere_cheating}")
        self.outputScaffolds(component_tuned_paths)
        return res
    
    #returns dict {old_id: new_str}
    def getAdditionalRephase(self, scaffolded_rephased, used_ids, paths):
        res = {}
        label_switch = {"HAPLOTYPE1": "HAPLOTYPE2", "HAPLOTYPE2": "HAPLOTYPE1"}
        prefix_switch = {"pat": "mat", "mat": "pat", "haplotype1": "haplotype2", "haplotype2": "haplotype1"}
        for nor_old_id in scaffolded_rephased:
            if not self.isHaploidPath(paths, nor_old_id):
                #very unoptimal but this should not be frequent
                homologous_paths = self.match_graph.getHomologousPaths(paths, nor_old_id)
                for alt_path_id in paths.getPathIds():
                    #scaffolded path should not be significantly shorter than homologous alternative to recolor
                    if paths.getLength(nor_old_id) * 2 > paths.getLength(alt_path_id) and alt_path_id in homologous_paths:
                    #and self.match_graph.isHomologousPath([paths.getPathById(nor_old_id), paths.getPathById(alt_path_id)], [paths.getLength(nor_old_id), paths.getLength(alt_path_id)]):
                        #alt_path_id is homologous to something rephased and thus should be rephased too.
                        old_label = paths.getLabel(alt_path_id)
                        splitted_id = alt_path_id.split("_")
                        old_prefix = splitted_id[0]
                        if alt_path_id in used_ids:
                            self.logger.info(f"Not rephasing {alt_path_id} because it is scaffolded anyway")
                            continue
                        if old_label in label_switch and old_prefix in prefix_switch:
                            new_label = label_switch[old_label]
                            new_prefix = prefix_switch[old_prefix]                            
                            new_id = new_prefix + "_" + "_".join(splitted_id[1:])
                            self.logger.warning(f"Rephasing {alt_path_id} to {new_id} because of {nor_old_id} scaffolded")
                            new_tsv = "\t".join([new_id, paths.getPathString(alt_path_id), new_label])
                            res[alt_path_id] = new_tsv
                        else:
                            self.logger.warning(f"Rephasing failed for {alt_path_id} initiated by {nor_old_id}")                            
        return res
    #if only two paths with one telomere missing remains in connected component and all other long nodes are in T2T and some conditions on distance and length then we connect
    def completeConnectedComponent(self, prefinal_paths):
        res = path_storage.PathStorage(self.upd_G)
  
        #stats for logging
        added_t2t = 0
        total_jumps = 0

        scaffolded_paths = set()
        components = list(nx.weakly_connected_components(self.upd_G))
        node_comp = {}
        path_comp = {}
        node_path = {}
  
        nor_components = list()
        used_nodes = set()
        for comp_id in range(0, len(components)):
            current = set()
            for node in components[comp_id]:
                nor_node = gf.nor_node(node)
                if nor_node in used_nodes:
                    continue
                else:
                    current.add(nor_node)
                    used_nodes.add(nor_node)
            if len(current) > 0:
                nor_components.append(current)
  
        for comp_id in range(0, len(nor_components)):
            for node in nor_components[comp_id]:
                node_comp[node] = comp_id
  

        self.logger.info(f"Checking {len(nor_components)} weakly connected components ")
        for path_id in prefinal_paths.getPathIds():
            cur_color = 0
            for node in prefinal_paths.getPathById(path_id):
                if node in self.upd_G.nodes:
                    nor_node = gf.nor_node(node)
                    if cur_color != 0 and node_comp[nor_node] != cur_color:
                        cur_color = -1
                    #some nodes in multiple paths but who cares
                    node_path[nor_node] = path_id                      
                    cur_color = node_comp[nor_node]
            path_comp[path_id] = cur_color

        for comp_id in range(0, len(nor_components)):
            valid = True
            comp_paths_ids = set()
            for node in nor_components[comp_id]:
                if node in node_path:
                    #if there is a path that jumps from component somewhere else, we do not want to scaffold anything
                    if node_path == -1:
                        valid = False
                    comp_paths_ids.add(node_path[node])
            if not valid:
                continue

            missing_telos = 0
            to_scaffold = []
            total_long_paths = 0
            telos = {}
            t2t_paths = 0

            for path_id in comp_paths_ids:
                if prefinal_paths.getLength(path_id) > ScaffoldGraph.MIN_EXPECTED_T2T:
                    total_long_paths += 1
                end_telo = False
                for n in self.upd_G.successors(prefinal_paths.getPathById(path_id)[-1]):
                    if n in self.tel_nodes:
                        end_telo = True                    
                start_telo = False
                for n in self.upd_G.predecessors(prefinal_paths.getPathById(path_id)[0]):
                    if n in self.tel_nodes:
                        start_telo = True
                telos[path_id] = (start_telo, end_telo)
                #component is not finished yet
                #path with only one telomere - we are interested
                if end_telo != start_telo:
                    to_scaffold.append(path_id)
                    missing_telos += 1
                else:
                    #Short non-telo paths are ignored, short t2t is weird case of telomeric-only node and also should be ignored
                    if prefinal_paths.getLength(path_id) > ScaffoldGraph.MIN_EXPECTED_T2T:
                        if start_telo and end_telo:
                            t2t_paths += 1
                        #component not finished, stopping
                        elif (not start_telo) and not (end_telo):
                            self.logger.debug(f"Found path {path_id} of length {prefinal_paths.getLength(path_id)} without telomeres, do nothing with component ")
                            missing_telos = 1000
                            break
                #only two non-t2t paths of reasonable length in connected component
            if (total_long_paths > 0):
                self.logger.debug(f"Component {comp_id} has {missing_telos} non-t2t telos, {t2t_paths} t2t paths, {total_long_paths} long paths")

            if missing_telos == 2:
                or_ids = []
                for i in range(0, 2):
                    if telos[to_scaffold[i]][i]:
                        or_ids.append(to_scaffold[i] + "+")
                    else:
                        or_ids.append(to_scaffold[i] + "-")

                dist = self.orPathIdDist(or_ids[0], or_ids[1], prefinal_paths, True)
                #Paths are not far away
                self.logger.info (f"Potential t2t connectivity cheat {or_ids} dist {dist}")    

                if dist < ScaffoldGraph.CLOSE_IN_GRAPH:
                    sum_len = prefinal_paths.getLength(to_scaffold[0]) + prefinal_paths.getLength(to_scaffold[1])                        
                    similar_len_t2t = 0
                    distant_len_t2t = 0
                    haploid_paths = 0
                    for alt_path_id in comp_paths_ids:
                        cur_len = prefinal_paths.getLength(alt_path_id) 
                        if self.isHaploidPath(prefinal_paths, alt_path_id):
                            haploid_paths += 1
                        if telos[alt_path_id][0] and telos[alt_path_id][1] and cur_len > ScaffoldGraph.MIN_EXPECTED_T2T:
                            if cur_len > sum_len * ScaffoldGraph.SIMILAR_LEN_FRAC and cur_len * ScaffoldGraph.SIMILAR_LEN_FRAC < sum_len:
                                similar_len_t2t += 1
                            else:
                                distant_len_t2t += 1              
                    #Either haploid or there is homologous T2T of similar length, not too short              
                    if ((distant_len_t2t == 0 and haploid_paths ==2) or similar_len_t2t > 0) and sum_len > ScaffoldGraph.MIN_EXPECTED_T2T:
                        node_sec = []                            
                        self.logger.info(f"t2t connectivity cheat worked {or_ids}!")
                        new_tsv_tuple = self.scaffoldFromOrPathIds(or_ids, prefinal_paths)
                        res.addPath(new_tsv_tuple[0])
                        if new_tsv_tuple[1]:
                            added_t2t += 1
                        total_jumps += new_tsv_tuple[2]
                        for i in range(0, 2):
                            scaffolded_paths.add(to_scaffold[i].strip('-+'))
                    else:
                        self.logger.info(f"t2t connectivity cheat {or_ids} failed, total_len {sum_len} distant len count {distant_len_t2t} haploid paths {haploid_paths} similar len {similar_len_t2t}")
                else:
                    self.logger.info(f"t2t connectivity cheat {or_ids} failed, too far {dist}")

        for path_id in prefinal_paths.getPathIds():
            if not (path_id in scaffolded_paths):
                res.addPath(prefinal_paths.getPathTsv(path_id))
        self.logger.warning (f"Total last telomere tuned scaffolds (all t2t) {added_t2t}, jumps {total_jumps}")
        return res                            
    
    #Returns (new_tsv_line, is_t2t, num_jumps, nor_rephased_ids)
    def scaffoldFromOrPathIds(self, or_ids, prefinal_paths):
        scf_path = []
        largest_label = "NA"
        largest_id = gf.nor_path_id(or_ids[0])
        largest_len = 0
        for i in range (0, len(or_ids)):
            nor_path_id = gf.nor_path_id(or_ids[i])
            if i != 0:
                #something non-default
                if or_ids[i][-1] == "+":
                    after_gap = prefinal_paths.getPathById(nor_path_id)[0]
                else:
                    after_gap = gf.rc_path(prefinal_paths.getPathById(nor_path_id))[0]
                prev_nor_id = gf.nor_path_id(or_ids[i-1])
                if or_ids[i - 1][-1] == "+":
                    before_gap = prefinal_paths.getPathById(prev_nor_id)[-1]
                else:
                    before_gap = gf.rc_path(prefinal_paths.getPathById(prev_nor_id))[-1]                
                gap_len = min(ScaffoldGraph.DEFAULT_GAP_SIZE, round(self.orNodeDist(before_gap, after_gap)))
                if gap_len != 0:
                    scf_path.append(f"[N{gap_len}N:scaffold]")
                else:
                    self.logger.info(f"Zero gap between {or_ids[i-1]} and {or_ids[i]}")
            if or_ids[i][-1] == "+":
                scf_path.extend(prefinal_paths.getPathById(nor_path_id))
            else:
                scf_path.extend(gf.rc_path(prefinal_paths.getPathById(nor_path_id)))
            #name and haplo of scaffolded path is taken from longest contig in scaffold
            if prefinal_paths.getLength(nor_path_id) > largest_len:
                largest_len = prefinal_paths.getLength(nor_path_id)
                largest_id = nor_path_id
                largest_label = prefinal_paths.getLabel(nor_path_id)
        rephased = []
        for or_id in or_ids:
            nor_path_id = gf.nor_path_id(or_id)
            label = prefinal_paths.getLabel(nor_path_id)
            if label != "NA" and label != largest_label:
                rephased.append(nor_path_id)
        #dangerous_nodes can be not correct on non-first iteration, but anyway debug only
        for i in range (0, len(or_ids) - 1):
            swap_info = (or_ids[i].strip('-+'), or_ids[i+1].strip('-+'), or_ids[i][-1] + or_ids[i+1][-1])
            if (swap_info) in self.dangerous_swaps:
                self.logger.warning (f"consecutive pair, {swap_info} did signficant changes on hi-c counts based on {self.dangerous_swaps[swap_info]}")
        path_str = "\t".join([largest_id, ",".join(scf_path), largest_label])
        telo_end = False
        telo_start = False
        for tel_node in self.tel_nodes:
            if self.upd_G.has_edge(scf_path[-1], tel_node):
                telo_end = True
            if self.upd_G.has_edge(tel_node, scf_path[0]):
                telo_start = True
        if len(or_ids) > 1:
            self.logger.info (f"Added SCAFFOLD {or_ids} {telo_start} {telo_end} ") 
        return (path_str, (telo_start and telo_end), len(or_ids) - 1, rephased)

#large haploid paths should be assigned based on connectivity. Complete human chrX will always be in hap1 
#TODO revisit after rukki's update!
    def fixHaploidNames(self, paths: path_storage.PathStorage):
        large_haploid_ids = []
        haploids = self.getHaploidPaths(paths)
        for path_id in haploids:
            if paths.getLength(path_id) > ScaffoldGraph.LONG_HAPLOID_CUTOFF:
                close_diploid_path = False
                for alt_path in paths.getPathIds():
                    if not (alt_path in haploids) and paths.getLength(alt_path) > ScaffoldGraph.LONG_HAPLOID_CUTOFF:
                        if self.norPathIdDist(alt_path, path_id, paths, False) < ScaffoldGraph.TOO_FAR:
                            close_diploid_path = True
                            break
                if not (close_diploid_path):
                    large_haploid_ids.append(path_id)
                else:
                    self.logger.warning(f"large haploid and diploid paths connected, {path_id} and {alt_path}, lens {paths.getLength(path_id)} {paths.getLength(alt_path)}, no relabeing")
        large_haploid_ids.sort(key=lambda x: paths.getLength(x), reverse=True)

        for path_id in paths.getPathIds():            
            if path_id.split("_")[0] == "mat":
                self.logger.info("Trio labeling detected, haploid haplotype reassignment cancelled")
                return
            
        names_prefix = {}
        names_prefix["HAPLOTYPE1"] = "haplotype1_from_"
        names_prefix["HAPLOTYPE2"] = "haplotype2_from_"


        large_haploid_info = []
        reassigned = 0
        reassigning_hap_id = 2
        not_connected_hap_assign = 0
        for i in range (0, len(large_haploid_ids)):
            path_id = large_haploid_ids[i]
            cur_label = paths.getLabel(path_id)
            #Two largest large haploids assigned to same haplotype are not normal            
            clothest_dist = ScaffoldGraph.TOO_FAR
            closest_id = -1
            for j in range (0, i):
                dist = self.norPathIdDist(large_haploid_ids[j], path_id, paths,False)
                if dist < clothest_dist:
                    clothest_dist = dist
                    closest_id = j
            #If we are close to some longer haploid use its precomputed label else alterate between two haplotypes
            if clothest_dist >= 2* ScaffoldGraph.CLOSE_IN_GRAPH: 
                new_label = "HAPLOTYPE" + str(reassigning_hap_id)
                reassigning_hap_id = 3 - reassigning_hap_id
                not_connected_hap_assign += 1
            else:
                new_label = large_haploid_info[closest_id][2]
            new_id = path_id
            #TODO should it be CLOSE_IN_GRAPH or just same connected component? Now with XY glued together by PAR we'll still have different colors                
            if new_label != cur_label:                                    
                utig_id = path_id.split("_")[-1]                    
                new_id = names_prefix[new_label] + utig_id
                self.logger.info(f"Reasigning haploid label for {path_id} {cur_label} to {new_id} {new_label}, dist {clothest_dist} len {paths.getLength(path_id)}")
                reassigned += 1                    
            else:                        
                self.logger.info(f"Reassignment of haploid label for {path_id} not required, dist {clothest_dist} len {paths.getLength(path_id)}")                
            large_haploid_info.append([new_id, ",".join(paths.getPathById(path_id)), new_label])            
        for path_id in large_haploid_ids:
            paths.removePath(path_id)
        for info in large_haploid_info:
            paths.addPath("\t".join(info))
        self.logger.info(f"Reassigned {reassigned} large haploids, groups of distant haploid paths detected: {not_connected_hap_assign}")
        if not_connected_hap_assign != 0 and not_connected_hap_assign != 2:
            self.logger.warning(f"Strange number of not connected haploid groups {not_connected_hap_assign}")


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
    #TODO: this works with hic_mapping.byread.output file that should be exterminated after phasing refactoring
    def get_connections_porec(self, alignment_file, use_multimappers:bool):
        res = {}
        unique_storage = {}
        #A01660:39:HNYM7DSX3:1:1101:1696:29982   utig4-73        utig4-1056      1       16949880        78591191
        ind = 0
        for line in open (alignment_file):
            ind += 1
            if (ind % 10000000 == 0):
                self.logger.debug (f"Processed {ind} pore-c alignments")

            arr = line.split()
            first = arr[1]
            second = arr[2]
            first_coords = int(arr[4])
            second_coords = int(arr[5])
            if first > second:
                first, second = second, first
                first_coords, second_coords = second_coords, first_coords
            weight = self.INT_NORMALIZATION
            names_pair = (first, second)
            coords_pair = (first_coords, second_coords)
            if not names_pair in res:
                res[names_pair] = []
            if not names_pair in unique_storage:
                unique_storage[names_pair] = []
            res[names_pair].append((coords_pair[0], coords_pair[1], weight))
            unique_storage[names_pair].append((coords_pair[0], coords_pair[1], weight))
        self.logger.debug (f"Finally processed {ind} pore-c alignments")
        return res, unique_storage
    
    def get_connections_bam(self, bam_filename, use_multimappers:bool):
        res = {}
        unique_storage = {}
        #A01660:39:HNYM7DSX3:1:1101:1696:29982   utig4-73        utig4-1056      1       16949880        78591191
        total_reads = 0
        unique_pairs = 0
        valid_pairs = 0
        all_pairs = 0
        
        created_pair = 0
        total_inserted = 0
        total_compressed = 0
        bamfile = pysam.AlignmentFile(bam_filename, "rb")
        cur_name = ""
        reads = []
        prev_read = None
        prev_name = ""
        
        for cur_read in bamfile:
            if (total_reads % 10000000 == 0):
                if (total_reads % 100000000 == 0):
                    before_compressed = total_compressed
                    self.logger.debug("Starting map compression...")
                    for names_pair in res:
                        if len(res[names_pair]) > 1000:
                            names_sorted = sorted(res[names_pair])
                            compressed = []
                            slen = len(names_sorted)
                            ii = 0
                            while ii < slen:
                                jj = ii + 1
                                cur_w = names_sorted[ii][2]
                                while jj < slen and names_sorted[jj][0] == names_sorted[ii][0] and names_sorted[jj][1] == names_sorted[ii][1]:
                                    cur_w += names_sorted[jj][2]
                                    jj += 1
                                    total_compressed += 1
                                compressed.append((names_sorted[ii][0], names_sorted[ii][1], cur_w))
                                ii = jj
                            res[names_pair] = compressed
                    self.logger.debug(f"On current iteration compressed {total_compressed - before_compressed}")
                self.logger.debug (f"Processed {total_reads} alignment strings")
                self.logger.debug (f"Of them unique vaild unique pairs {unique_pairs}, total pairs {all_pairs} total valid {valid_pairs} ")
                self.logger.debug (f"Current memory usage {(psutil.virtual_memory().used / 1024 / 1024 / 1024):.2f} GB")
                self.logger.debug (f"Created new pairs {created_pair} total inserted pairs = {total_inserted} total compressed = {total_compressed}")
                #gc.collect()
                #self.logger.debug (f"Current memory usage after GC {psutil.virtual_memory().used / 1024 / 1024 / 1024}GB")
                #self.logger.debug (f"Mem usage of main map {asizeof.asizeof(res) / 1024 / 1024 / 1024}GB")

                #if total_reads == 20000000:
                #    exit()
            total_reads += 1
            cur_name = cur_read.query_name
            if cur_name == prev_name:
                #TODO: poreC is not compatible with this now
                reads = (prev_read, cur_read)
#                  if read.is_paired:
                all_pairs += 1
                if prev_read.mapping_quality > 0 and cur_read.mapping_quality > 0:
                    #TODO: special storage possibly not needed, just check weights?
                    unique_pairs += 1
                    #if not (prev_read.reference_name, cur_read.reference_name) in unique_storage:
                    #    unique_storage[(prev_read.reference_name, cur_read.reference_name) ] = 0
                    #TODO: possibly require node1 < node2?
                    #next = (prev_read.reference_start, cur_read.reference_start, self.INT_NORMALIZATION) 
                    #unique_storage[(prev_read.reference_name, cur_read.reference_name) ]+= self.INT_NORMALIZATION

                names = [[prev_read.reference_name], [cur_read.reference_name]]
                coords = [[prev_read.reference_start], [cur_read.reference_start]]
                quals = [prev_read.mapping_quality, cur_read.mapping_quality]
#                    names = read.reference_name
                valid = True
                if use_multimappers:
                    i = 0
                    for read in reads:
                        if read.has_tag("XA"):
                            nm = int(read.get_tag("NM"))
                            if read.has_tag("XS") and read.get_tag("XS") == read.get_tag("AS"):  
                                for xa in read.get_tag("XA")[:-1].split(";"):
                                    xa_arr = xa.split(",")
                                    #TODO: remove later
                                    for pname in names[1 - i]:
                                        if xa_arr[0] == pname:
                                            valid = False
                                    #do not want to do cigar comparsion
                                    if int(xa_arr[3]) == nm:
                                        names[i].append(xa_arr[0])
                                        coords[i].append(int(xa_arr[1][1:]))
                        #Too many alignments, not reported in XA
                        #TODO: likely is prefiltered, check
                        elif read.mapping_quality == 0:
                            valid = False
                        i += 1
                    weight = self.INT_NORMALIZATION  // (len(names[0]) * len(names[1]))  

                    filtered_names = names
                    filtered_coords = coords
                    lname0 = len(filtered_names[0])
                    lname1 = len(filtered_names[1])
                    #self.logger.debug (f" {valid} {lname0} {lname1}")
                    if valid and lname0 > 0 and lname1 > 0:
                        valid_pairs += 1
                        for i in range (0, lname0):
                            node_f_len = self.uncompressed_lens[filtered_names[0][i]]
                            #TODO: remove check, this memory opt is already in hic_prefilter
                            if node_f_len < ScaffoldGraph.SHORT_INGORED_NODE or (filtered_coords[0][i] > ScaffoldGraph.NEAR_PATH_END and node_f_len - filtered_coords[0][i] > ScaffoldGraph.NEAR_PATH_END):
                                continue
                            for j in range (0, lname1):
                                node_s_len = self.uncompressed_lens[filtered_names[1][j]]
                                if node_s_len < ScaffoldGraph.SHORT_INGORED_NODE or (filtered_coords[1][j] > ScaffoldGraph.NEAR_PATH_END and node_s_len - filtered_coords[1][j] > ScaffoldGraph.NEAR_PATH_END):
                                    continue
                                if (filtered_names[0][i] < filtered_names[1][j]):
                                    names_pair = (filtered_names[0][i], filtered_names[1][j])
                                    pre_coords_pair = (filtered_coords[0][i], filtered_coords[1][j])
                                else :
                                    names_pair = (filtered_names[1][j], filtered_names[0][i])
                                    pre_coords_pair = (filtered_coords[1][j], filtered_coords[0][i])
                                #To round to closest integer and not always down
                                coords_pair = tuple(((c + self.APPROXIMATE_COORDS_HALF) // self.APPROXIMATE_COORDS) * self.APPROXIMATE_COORDS for c in pre_coords_pair)
                                if not names_pair in res:
                                    res[names_pair] = []
                                    created_pair += 1
                                res[names_pair].append((coords_pair[0], coords_pair[1], weight))                            
                                total_inserted += 1
                                
                                if quals[0] > 0 and quals[1] > 0:
                                    if not names_pair in unique_storage:
                                        unique_storage[names_pair] = []
                                    unique_storage[names_pair].append((coords_pair[0], coords_pair[1], weight))                            
                                #TODO: possibly require node1 < node2?
                                #next = (coords[0][i], coords[1][j], weight)
                                #res[(names[0][i], names[1][j])].append(next)

                                #no more check on homologous regions?
                                #res[(names[0][i], names[1][j])] += weight
            prev_read = cur_read
            prev_name = cur_name
        return res, unique_storage

    #TODO: move to matchGraph
    
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
                #with multimappers we likely dooo...
                if self.rukki_paths.getLength(path_ids[i]) <= self.TOO_FAR:
                    correct_or = ""
                    if self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1]:
                        correct_or = correct_orientations[i]
                    elif self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][1] and not self.scaffold_graph.nodes[path_ids[i]+"+"]['telomere'][0]:
                        correct_or = gf.rc_orientation(correct_orientations[i])
                    if correct_or != "":            
                        pair_or = ["",""]
                        pair_or[i] = correct_or                        
                        self.logger.debug (f"tuning pair {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                        for second_or in ('-', '+'):
                            correct_pair = pair_or.copy()                       
                            correct_pair[1 - i] = second_or
                            incorrect_pair = correct_pair.copy()
                            incorrect_pair[i] = gf.rc_orientation(correct_or)
                            correct_pair = "".join(correct_pair)
                            incorrect_pair = "".join(incorrect_pair)                            
                            #just logging
                            if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS* self.INT_NORMALIZATION and scores[incorrect_pair]  > scores[correct_pair]:
                                self.dangerous_swaps[(path_ids[0], path_ids[1], correct_pair)] = "telomeres"
                                self.logger.debug(f"Dangerous telomeric tuning pair {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']} from {incorrect_pair} to {correct_pair}")                    
                            self.logger.debug (f"moving {incorrect_pair} to {correct_pair}")
                            if self.rukki_paths.getLength(path_ids[i]) <= self.NEAR_PATH_END:
                                scores[correct_pair] += scores[incorrect_pair]
                            else:
                                if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS * self.INT_NORMALIZATION and scores[incorrect_pair]  > scores[correct_pair]:
                                    self.logger.debug(f"Dangerous telomeric tuning pair too long to move {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                            scores[incorrect_pair] = 0                            
                        self.logger.debug (f"telomeric tuned pair {path_ids}, scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}") 
        
            for i in range (0, 2):
                #TODO: WIP
                for fixed_orientation in ('-', '+'):
                    shortest_paths = {'-':1000000000, '+':1000000000}
                    min_cutoff = min(self.CLOSE_IN_GRAPH, self.rukki_paths.getLength(path_ids[i]) / 4)                        
                    max_cutoff = self.rukki_paths.getLength(path_ids[i]) * 3 / 4
                    for orient in ('-', '+'):
                        to_check_ids =["",""]
                        to_check_ids[1 - i] = path_ids[1 - i] + fixed_orientation
                        to_check_ids[i] = path_ids[i] + orient
                        shortest_paths[orient] = self.orPathIdDist(to_check_ids[0], to_check_ids[1], self.rukki_paths, True)
                        self.logger.debug(f"Checking dists {to_check_ids} index {i} dist {shortest_paths[orient]} cutoffs {min_cutoff} {max_cutoff}")

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
                        incorrect_pair[i] = gf.rc_orientation(correct_or)
                        incorrect_pair[1 - i] = fixed_orientation   
                        correct_pair = "".join(correct_pair)
                        incorrect_pair = "".join(incorrect_pair) 
                        #just logging
                        if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS * self.INT_NORMALIZATION and scores[incorrect_pair] > scores[correct_pair]:
                            self.dangerous_swaps[(path_ids[0], path_ids[1], correct_pair)] = "connectivity"
                            self.logger.debug(f"Dangerous connectivity tuning pair {path_ids}, i {i} scores {scores}from {incorrect_pair} to {correct_pair}")                    
                        if self.rukki_paths.getLength(path_ids[i]) <= self.NEAR_PATH_END:
                            scores[correct_pair] += scores[incorrect_pair]
                        else:
                            if scores[incorrect_pair] >= ScaffoldGraph.MIN_LINKS * self.INT_NORMALIZATION and scores[incorrect_pair]  > scores[correct_pair]:
                                self.logger.debug(f"Dangerous connectivity  tuning pair too long to move {path_ids}, i {i} scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}")                    
                        #TODO: this may happen twice or once!!!
                        #scores[correct_pair] *= self.CONNECTIVITY_MULTIPLICATIVE_BONUS
                        scores[incorrect_pair] = 0 
                        #self.logger.debug (f"Connectivity tuned pair {path_ids}, scores {scores}, tels {self.scaffold_graph.nodes[path_ids[i]+'+']['telomere']}") 
            #Code duplication:(
            for first_or in ('-', '+'):
                for second_or in ('-', '+'):
                    or_str = first_or + second_or
                    to_check = [path_ids[0]+first_or, path_ids[1]+second_or]
                    dist_paths = self.orPathIdDist(to_check[0], to_check[1], self.rukki_paths, True)     
                    if dist_paths < self.CLOSE_IN_GRAPH and dist_paths < self.rukki_paths.getLength(path_ids[0]) / 4 and dist_paths < self.rukki_paths.getLength(path_ids[1]) / 4:
                        scores[or_str] *= self.CONNECTIVITY_MULTIPLICATIVE_BONUS
                        self.logger.debug(f"Connectivity bonus applied pair {to_check_ids}, scores {scores}")

        return scores

    #return scores for each of the orientations, ++, -+, +-, --,
    #orientation within the path, it is not changing!
    def getNodePairConnections(self, pair, orientations, connections, shift_before, shift_after, lens):    
        #This is scores for PATH orientation
        scores = {"++":0, "-+":0, "+-":0, "--":0, "middle":0}    
        filtered = 0
        not_filtered = 0
        if self.match_graph.hasEdge(pair[0], pair[1]):
            intervals = self.match_graph.getEdgeAttribute(pair[0], pair[1], 'intervals')
        else:
            intervals = [[],[]]
        #Force ignoring links between homologous nodes
        if self.match_graph.isHomologousNodes(pair[0], pair[1], True):
            return scores
        
        #res[(names[1][j], names[0][i])].append([coords[1][j], coords[0][i], weight]) 
        rc_pair = (pair[1], pair[0])

        cons = []

        #stored only half
        #TODO: not optimal
        if rc_pair in connections:
            for conn in connections[rc_pair]:
                cons.append((conn[1], conn[0], conn[2]))
        if pair in connections:
            for conn in connections[pair]:        
                cons.append(conn)

        for conn in cons:        
            in_homo = False
            for i in range (0, len(intervals[0])):
                local_homo = True
                for j in range (0, 2):
                    #coords were approximated so adjusting. not sure we'll actually need it
                    if conn[j] <  intervals[j][i][0] + self.APPROXIMATE_COORDS_HALF or conn[j] > intervals[j][i][1] - self.APPROXIMATE_COORDS_HALF:
                        local_homo = False
                    #TODO: check!!
                    if local_homo:
                        in_homo = True
                        break
            if in_homo:
                filtered += 1
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
                scores[str] += conn[2]
            else:
                scores["middle"] += conn[2]
        self.logger.debug (f"Scores for {pair} {scores} {orientations} filtered/not_filtered {filtered} {not_filtered}")
        #self.logger.debug (f"Shifts {shift_before} {shift_after}")
        #self.logger.debug (f"Ignored intervals {intervals}")
        return scores
    
    #TODO: remove, not used anymore
    def getInterestingNodes(self, lens):
        interesting = set()
        for path_id in self.rukki_paths.getPathIds():            
            #TODO possibly reconsider this condition; /2 because of compressed/uncompressed
            if self.rukki_paths.getLength(path_id) < ScaffoldGraph.MIN_PATH_TO_SCAFFOLD /2 :
                continue
            total_len = 0
            for or_node in self.rukki_paths.getPathById(path_id):
                nor_node = or_node.strip('-+')
                if nor_node in lens and lens[nor_node] > ScaffoldGraph.SHORT_INGORED_NODE and or_node in self.multiplicities and self.multiplicities[or_node] == 1:
                    total_len += lens[nor_node]
            before = 0
            after = total_len
            for or_node in self.rukki_paths.getPathById(path_id):
                nor_node = or_node.strip('-+')
                if nor_node in lens and lens[nor_node] > ScaffoldGraph.SHORT_INGORED_NODE and or_node in self.multiplicities and self.multiplicities[or_node] == 1:
                    after -= lens[nor_node]
                    if before < ScaffoldGraph.NEAR_PATH_END or after < ScaffoldGraph.NEAR_PATH_END:
                        interesting.add(int(nor_node[6:]))
                    before += lens[nor_node]
        return frozenset(interesting)
    
    #Just to reduce debug flood
    def isPathPairForDebug(self, path_ids, path_storage):
        if path_storage.getLength(path_ids[0]) < ScaffoldGraph.MIN_PATH_TO_SCAFFOLD and path_storage.getLength(path_ids[1]) < ScaffoldGraph.MIN_PATH_TO_SCAFFOLD:
            return False
        else:
            return True
        
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
                #TODO: should we really ignore short nodes?
                if not (nor_s in lens) or lens[nor_s] < ScaffoldGraph.SHORT_INGORED_NODE:
                    continue
                if not second in self.multiplicities or self.multiplicities[second] > 1:
                    continue
                pair = (nor_f, nor_s)
                orientations = (first[-1], second[-1])
                if not pair in connections and not (pair[1], pair[0]) in connections:
                    continue
                shift_before = []
                shift_after = []
                for i in range (0, 2):
                    shift_before.append(length_before[i][pair[i]])
                    shift_after.append(total_len[i] - shift_before[i] - lens[pair[i]])
                scores_pair = self.getNodePairConnections(pair, orientations, connections, shift_before, shift_after, lens)
                
                for key in scores_pair:
                    scores[key] += scores_pair[key] 
        if self.isPathPairForDebug(path_ids, self.rukki_paths):
            self.logger.debug (f"Pathscores for {path_ids} {scores}")
        scores = self.fixOrientation(path_ids, scores)
        for orientation in scores:
            if len(orientation) == 2:
                if self.isNextByRef(path_ids[0] + orientation[0], path_ids[1] + orientation[1]):
                    self.logger.debug (f"Reference connection found! {path_ids} {orientation}")
                    scores[orientation] *= self.REFERENCE_MULTIPLICATIVE_BONUS
            scores[orientation] /= self.INT_NORMALIZATION

            very_short = 0
            very_long = 0
            for i in range (0, 2):
                if self.rukki_paths.getLength(path_ids[i]) > ScaffoldGraph.NEAR_PATH_END:
                    very_long += 1
                #constant tweeeeeing
                if self.rukki_paths.getLength(path_ids[i]) < ScaffoldGraph.NEAR_PATH_END/5:
                    very_short += 1
             #TODO: possibly do not need reordering from "wrong" end and compare with pathLength/2 after reordering removed?
            #Never want to apply this bonus twice
            if very_long <= 1 and very_short == 0:
                coeff = self.NEAR_PATH_END / min(self.rukki_paths.getLength(path_ids[0]), self.rukki_paths.getLength(path_ids[1]))
                self.logger.debug(f"Normalization is applied to {path_ids} , lens {self.rukki_paths.getLength(path_ids[0])} {self.rukki_paths.getLength(path_ids[1])} coefficent {coeff}")  
                scores[orientation] *= coeff
        return scores

#Do we need it for reruns?    
    def loadPresavedScores(self, presaved_file):
        pass

    #returns dict, {id:[present_start_relo, present_end_telo]}
    def getTelomericEnds(self):
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
            res[id] = [tel_start, tel_end]                
        return res

    #For each paths returns its connected component. If multiple, then 0
    def getPathColors(self, rukki_paths, graph):
        components = sorted(nx.weakly_connected_components(self.upd_G),key=len, reverse=True)
        node_colors = {}
        path_colors = {}
        for i in range (0, len(components)):
            for node in components[i]:
                node_colors[node] = i + 1
        for path_id in rukki_paths.getPathIds():
            current_colors = set()
            for node in rukki_paths.getPathById(path_id):
                current_colors.add(node_colors[gf.nor_node(node)])
            if len (current_colors) == 1:
                path_colors[path_id] = current_colors.pop()
            else:
                path_colors[path_id] = 0
        return path_colors
    
    
