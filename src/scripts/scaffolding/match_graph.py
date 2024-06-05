#!/usr/bin/env python3
import networkx as nx
import graph_functions as gf
import sys
from scaffolding import logger_wrap

#classes for matchGraph construction
class HomologyInfo:
    def __init__(self, node1, node2, len1, len2):
        self.nodes = [node1, node2]
        self.len = [len1, len2]
        self.covered = [0, 0]
        #two lists for each nodes, each list constains [st_i, 1] and [end_i, -11]
        self.intervals = [[],[]]

        #orientation defined as orientation of the largest match
        self.largest_interval = 0
        self.largest_interval_center = [0, 0]
        self.orientation = ""

        #TODO: should be used for homology check in scaffolding, with specific IDY cutoff and not sorted
        self.filtered_intervals = [[],[]]        

    def parseIDY(self, idy):
        return float(idy.split(":")[2])
    
    def addInterval(self, intervals, orientation, idy):
        real_idy = self.parseIDY(idy)
        #Constant not depending on genome, intervals too similar for hi-c alignment to use
        if real_idy > 0.995:
            for i in range(0, 2):
                self.filtered_intervals[i].append(intervals[i])

        #do we need to check here bad matches?
        for i in range(0, 2):
            self.intervals[i].append([intervals[i][0],1])
            self.intervals[i].append([intervals[i][1],-1])
                
        int_len = min(intervals[0][1] - intervals[0][0], intervals[1][1] - intervals[1][0])
        if self.largest_interval < int_len:
            self.largest_interval = int_len
            self.orientation = orientation
            self.largest_interval_center = [(intervals[0][1] + intervals[0][0])/2, (intervals[1][1] + intervals[1][0])/2]

    def fillCoverage(self):
        for i in range(0, 2):
            total_c = 0
            self.intervals[i].sort()
            state = 0
            prev = 0
            for coord_pair in self.intervals[i]:
                if state > 0:
                    total_c += coord_pair[0] - prev
                prev = coord_pair[0]
                state += coord_pair[1]
            self.covered[i] = total_c   
    
    def getCoveredLen(self):
#in weird case matches can be larger than seq, avoiding
        return min(self.covered[0], self.covered[1], self.len[0], self.len[1])
    
    def getMinLength(self):
        return min(self.len[0], self.len[1])

class HomologyStorage:
    #{node1: {node2: HomologyInfo(node_1, node_2)}}
    def __init__(self, logger, mashmap_file, min_alignment):
        self.homologies = {}   
        self.lens = {}

        self.logger = logger_wrap.UpdatedAdapter(logger, self.__class__.__name__)

        total_lines = 0
        used_lines = 0
        for line in open(mashmap_file, 'r'):
            arr = line.strip().split()
            total_lines += 1
            if len(arr) < 11:
                continue
            if int(arr[10]) < min_alignment:
                continue
            #self mapping
            if arr[0] == arr[5]:
                continue
            if len(arr[0]) < 3:
                continue
            used_lines +=1
            #utig4-0 2145330 0       990000  +       utig4-0 2145330 12      994065  37      994053  51      id:f:0.999992   kc:f:0.874893
            self.addHomology(arr[0], arr[5], int(arr[1]), int(arr[6]), [[int(arr[2]), int(arr[3])], [int(arr[7]), int(arr[8])]], arr[4], arr[12])
        self.logger.info(f"Loaded {used_lines} out of {total_lines} mashmap lines")
        self.logger.info(f"{len(self.homologies)} nodes have at least one used homology")
        self.fillCoverage()

#do we want to add other direction of pair?
    def addHomology(self, node1, node2, len1, len2, intervals, orientation, idy):
        if not node1 in self.homologies:
            self.homologies[node1] = {}
        if not node2 in self.homologies[node1]:
            self.homologies[node1][node2] = HomologyInfo(node1, node2, len1, len2)
        self.homologies[node1][node2].addInterval(intervals, orientation, idy)
        self.lens[node1] = len1
        self.lens[node2] = len2
    #covered length of homology
    def fillCoverage(self):
        for node1 in self.homologies:
            for node2 in self.homologies[node1]:    
                self.homologies[node1][node2].fillCoverage()

    def isValid(self, node1, node2):
        if node1 in self.homologies and node2 in self.homologies[node1]:
            return True
        else:
            return False
    #extracting length, sometimes we do not haev this info in other places
    def getLength(self, node):
        return self.lens[node]
        
# homology_weight - large negative something, min_big_homology - minimal length of homology to be taken in account. 
# Some shorter ones can still sometimes be used if they are in regular bulge_like structure
# G can be both directed and undirected graph
class MatchGraph:
    def __init__(self, mashmap_sim, G, homology_weight, min_big_homology, min_alignment, logger):
        self.matchGraph = nx.Graph()
        self.hom_storage = HomologyStorage(logger, mashmap_sim, min_alignment)
        self.G = G
        self.logger = logger_wrap.UpdatedAdapter(logger, self.__class__.__name__)

        
        #Do not want to think whether we use diGraph or Graph
        neighbours = {}
        indirect_nodes = set()
        for node in G.nodes():
            neighbours[node.strip('-+')] = set()
            indirect_nodes.add(node.strip('-+'))
        for edge in G.edges():
            for i in range(0, 2):
                neighbours[edge[i].strip('-+')].add(edge[1 - i].strip('-+'))
        self.logger.info(f"Loaded {len(self.hom_storage.homologies)} homologies")            
        for node1 in self.hom_storage.homologies:
            for node2 in self.hom_storage.homologies[node1]:
                #we deleted some nodes after mashmap
                #sys.stderr.write(f"Checking {node1} {node2} {self.hom_storage.homologies[node1][node2].getMinCovered()} {min_big_homology}\n")
                if node1 in indirect_nodes and node2 in indirect_nodes:
                    cur_homology = self.hom_storage.homologies[node1][node2].getCoveredLen()
                    if cur_homology > min_big_homology:
                        self.logger.debug(f"Adding normal edge {node1} {node2} {cur_homology}")
                        self.matchGraph.add_edge(node1, node2, homology_len = cur_homology)
                    else:
                        #less strict condition for bulge-like structure
                        #common = neighbours[0].intersection(neighbours[1])
                        len_cutoff = min(self.hom_storage.homologies[node1][node2].getMinLength(), min_big_homology)/2
                        # possibly check whether we have something already phased nearby?
                        # strict bulge-like condition, either R1[A/B]R2 or R1[A/B] or [A/B]R2. Digraph can be more precise here but anyway
                        if cur_homology > len_cutoff and neighbours[node1] == neighbours[node2] and len(neighbours[node1]) > 0:
                            self.logger.debug(f"Adding bulge-like edge {node1} {node2} {cur_homology} {len_cutoff}")
                            self.matchGraph.add_edge(node1, node2, homology_len = cur_homology)
    
            # while we build the initial partition give a big bonus edge for putting the homologous nodes into different partitions             
            # Adding an edge that already exists updates the edge data (in networkX graph)
            # At least one in the pair is the best similarity match for other (and second best is not too close to the best)
            # Consecutive edges are used to check but never assigned to be homologous
        for ec in self.matchGraph.edges():
            clear_best_match = False
            #homology storage may be asymetric, mashmap do not guararntee anything        
            #possibly should forcely symmetrize... 
            if self.hom_storage.isValid(ec[0], ec[1]) and (not G.has_edge(ec[0], ec[1])):
                for i in range (0, 2):
                    best_homology = True
                    best_len = self.matchGraph.edges[ec]['homology_len']
                    for adj_node in self.matchGraph.neighbors(ec[i]):
                        if adj_node != ec[1 - i] and  best_len * 0.8 < self.matchGraph.edges[ec[i], adj_node]['homology_len']:
                            best_homology = False
                    clear_best_match = clear_best_match or best_homology
                if clear_best_match:
                    self.matchGraph.add_edge(ec[0], ec[1], weight = homology_weight, homology_len = best_len, intervals = self.hom_storage.homologies[ec[0]][ec[1]].filtered_intervals, 
                                        orientation = self.hom_storage.homologies[ec[0]][ec[1]].orientation)
                else:
            #not really look like homologous node pair but still suspicious, lets just wipe the hi-c links but not prioritize splitting them to different partitions.                                
                    self.matchGraph.add_edge(ec[0], ec[1], weight = 0, homology_len = best_len, intervals = self.hom_storage.homologies[ec[0]][ec[1]].filtered_intervals,
                                        orientation = self.hom_storage.homologies[ec[0]][ec[1]].orientation)
            else:
                self.matchGraph.remove_edge(ec[0],ec[1])

        self.logger.info("Loaded match info with %d nodes and %d edges\n" % (self.matchGraph.number_of_nodes(), self.matchGraph.number_of_edges()))

        for d in sorted(self.matchGraph.edges()):
            self.logger.debug(f"homology edge {d} : {self.matchGraph.edges[d]}")

    def getMatchGraph(self):
        return self.matchGraph

#TODO: duplication?
    def isDiploid(self, component):
        hom_length = 0
        total_length = 0
        for node in component:
            total_length += self.G.nodes[node]['length']
            if node in self.matchGraph.nodes():
                for adj_node in self.matchGraph.neighbors(node):
                    if adj_node in component and self.matchGraph.edges[node, adj_node]['weight'] < 0:
                        hom_length += self.matchGraph.edges[node, adj_node]['homology_len']
        if hom_length * 2 > total_length:
            return True
        else:
            return False
        

    def isHomologous (self, paths, lens):
        hom_size = 0
        for p0 in paths[0]:
            for p1 in paths[1]:
                nor_p0 = p0.strip("-+")
                nor_p1 = p1.strip("-+")
                if self.matchGraph.has_edge(nor_p0, nor_p1):
                    if self.matchGraph.edges[nor_p0, nor_p1]['weight'] < 0:
                        
                        hom_size += self.matchGraph.edges[nor_p0, nor_p1]['homology_len']
        if hom_size * 2> lens[0]  or hom_size * 2> lens[1]:
            self.logger.debug (f"Found homologous paths {paths} with homology size {hom_size}")
            return True
        else:
            return False