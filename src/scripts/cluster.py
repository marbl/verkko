import logging
import sys
import random
import networkx as nx
import math
import os
from networkx.algorithms import community
import graph_functions as gf
from scaffolding import match_graph, logger_wrap


#TODO: move constants to some more relevant place
MIN_LEN = 200000  # best result so far with 200000

#likely do not need
#MAX_GRAPH_DIST = 100000000  # hi-c links over 10M are believed to be useless

KLIN_STARTS = 30  # number of different starts of kernighan lin
KLIN_ITER = 10000  # number of iterations inside kernighan lin

MAX_SHORT_COMPONENT = 100 # we remove from consideration each connected compoents of edges < MIN_SHORT_LEN that is larger than MAX_SHORT_COMPONENT
MAX_SHORT_LEN=200000 #we remove from consideration each connected compoents of edges < MIN_SHORT_LEN that is larger than
MIN_WEIGHT = 10 # Ignore edges with few links
SIGNIFICANT_MAJORITY = 2.5

RESULT_FILENAME = "hicverkko.colors.tsv"
FILTERED_GFA_FILENAME = "filtered.gfa"
HIC_COMPRESSED_FILENAME = "hic.byread.compressed"

CLEAR_HOMOLOGY = 250000 #for clear bulge like structures this limmit is decreased
MIN_ALIGNMENT = 50000 #smaller alingments will be filtered out
MAX_COV = 100  # tempora# ry coverage cutoff, currently replaced by median coverage from gfa
FIXED_WEIGHT = 100000  # best result so far with 100000 #currently replaced with max pairwise weight among datasets

MAX_RDNA_COMPONENT = 10000000 # maximal size of rDNA component, used for filtering out rDNA cluster only
MIN_RDNA_COMPONENT = 500000 
logger = logger_wrap.initLogger("phasing.log")


def check_non_empty(part, G):
    for p in part:
        if p in G.nodes:
            return True
    return False

#Currently not in use anymore
def IsTip(node, edges):
    for pref in ['>', '<']:
        ornode = pref + node
        if ornode not in edges or len(edges[ornode]) == 0:
            if gf.revnode(ornode) not in edges:
                continue
#rc to the tip predeccor
            for edge in edges[gf.revnode(ornode)]:
#alternative edge to tip that should be not a deadend
                for alt_tip in edges[gf.revnode(edge)]:
                    if alt_tip in edges and len(edges[alt_tip])> 0:
                        print (f'Tip {node}')
                        return True
    return False

#collapse a node, add links from
def collapseOrientedNode (edges, node):
    for pref in ['>', '<']:
        ornode = pref + node
        if not ornode in edges:
            continue
        if not gf.revnode(ornode) in edges:
            continue
        outgoing = set(edges[ornode])
        incoming = set()
        for n in edges[gf.revnode(ornode)]:
            incoming.add(gf.revnode(n))
        for inc in incoming:
            for outg in outgoing:
                edges[inc].add(outg)
    for pref in ['>', '<']:
        ornode = pref + node
        if not gf.revnode(ornode) in edges:
            continue
        incoming = set()
        for n in edges[gf.revnode(ornode)]:
            incoming.add(gf.revnode(n))
        for inc in incoming:
            edges[inc].discard(ornode)
    for pref in ['>', '<']:
        ornode = pref + node
        if ornode in edges:
            edges.pop(ornode)

#returns pair of nodes near PAR that should be connected via fake homology, currently placeholder
#TODO: this is likely subject for reconsidering when we'll move to other species

def checkXYcomponent(current_component, matchgraph, G, edges):
#minimal comp size. XY is vulnerable to gaps, so not 150M
    MIN_COMP_SIZE = 80000000
#maximal homozygous, PAR region
    MAX_HOM_LEN = 5000000     
#cutoff for additional fake homozygocity
    MIN_NEAR_PAR = 3000000
#TODO these constants are ABSOLUTELY human focused, and even if we'll save the function should be reconsidered.
#TODO: right way with new matchGraph...
    total_l = 0
    total_hom_len = 0
    for node in current_component:
        if node in G.nodes:
            total_l += G.nodes[node]['length']
            if node in matchgraph:
                total_hom_len += G.nodes[node]['length']

    if total_hom_len == 0 or total_hom_len > MAX_HOM_LEN or total_l < MIN_COMP_SIZE:
        return [0, 0]
    fake_hom = []
    for f_node in current_component:
        for s_node in current_component:
            if f_node < s_node and G.nodes[f_node]['length'] > MIN_NEAR_PAR and G.nodes[s_node]['length'] > MIN_NEAR_PAR:
                for fp in ['<', '>']:
                    for sp in ['<', '>']:
                        f_or = fp + f_node
                        s_or = sp + s_node
                        if not (f_or in edges) or not (s_or in edges):
                            continue
                        for nxt in edges[f_or]:
                            if nxt in edges[s_or]:
                                if len(fake_hom) == 0 or fake_hom[-1][0] != f_node or fake_hom[-1][1] != s_node:
                                    fake_hom.append([f_node, s_node])
    if len(fake_hom) > 1:
        logger.warning("MULTIPLE CANDIDATES FOR FAKE HOMOLOGY IN XY\n")
        for pair in fake_hom:
            logger.warning(f"{pair[0]} --- {pair[1]}\n")
        return [0, 0]
    elif len(fake_hom) == 1:
        logger.info(f"Adding FAKE HOMOLOGY between {fake_hom[0][0]} --- {fake_hom[0][1]}")
        return [fake_hom[0][0], fake_hom[0][1]]
    return [0, 0]
        
def fixUnbalanced(part, C, G):
    auxs = [0, 0]
    lens = [0, 0]
    for ind in range (0, 2):
        for node in part[ind]:
            if node.startswith('Aux'):
                auxs[ind] += 1
            elif node in G.nodes:
                lens[ind] += G.nodes[node]['length']
            else:
                print(f"Someting went wrong with node {node}")
#    print (f"In fix unbalanced {lens[0]} {lens[1]}")
    #just _a_ limitation on the swap length is required
    maxswap = min(lens[0]/10, lens[1]/10)
    curswap = 0
    edge_swapped =0
    for ind in range(0, 2):
        #do not move edges _to_ the part where aux edge is present -otherwise it could be swapped with aux
        if auxs[1 - ind] == 0:
#            print(f'processing components {part[0]} {part[1]}')
            changed = True
            while changed:
                changed = False
                for node in part[ind]:
                    if node in G:
                        cost = 0
                        for i in part[ind]:
                            if [i, node] in C.edges:
                                cost += C.edges[i, node]['weight']
                        for j in part[1 - ind]:
                            if [j, node] in C.edges:
                                cost -= C.edges[j, node]['weight']
                        if cost < 0 and curswap + G.nodes[node]['length'] < maxswap:
                            changed = True
                            curswap += G.nodes[node]['length']
                            part[1 - ind].add(node)
                            part[ind].remove(node)
                            edge_swapped +=1
                            print (f"SWAPPED {node}")
                            break
    if curswap > 0:
        print (f"Fixing uneven component, moved {curswap} bases and {edge_swapped} nodes")

def getMedianCov(nodeset):
    med_cov = -1
    total_length = 0
    for node in nodeset:
        total_length += node[1]['length']
    sum_l = 0
    sorted_nodes = sorted(nodeset, key=lambda node: node[1]['coverage'])
    for node in sorted_nodes:
        logger.debug(f'Node {node[0]} coverage {node[1]["coverage"]} length {node[1]["length"]}')
        sum_l += node[1]['length']
        if 2*sum_l > total_length:
            med_cov = node[1]['coverage']
            logger.info(f'Median coverage is {med_cov}')
            break
    return med_cov

def report_phased_node(node_id, hap, output_file):
    if hap == 1:
        output_file.write(f'{node_id}\t0\t100000\t0:100000\t#8888FF\n')
    elif hap == 0:
        output_file.write(f'{node_id}\t100000\t0\t100000:0\t#FF8888\n')

def create_initial_partition_noKL(C, matchGraph, seed):
    random.seed(seed)
    neg_graph = nx.Graph()
    for ec in matchGraph.edges():
        if ec[0] in C and ec[1] in C:
            if matchGraph.edges[ec]['weight'] < 0:
                for i in range (2):
                    if not ec[i] in neg_graph.nodes():
                        neg_graph.add_node(ec[i])
                neg_graph.add_edge(ec[0], ec[1])
        elif (ec[0] in C or ec[1] in C) and matchGraph.edges[ec]['weight'] < 0:
            logger.warning(f"Weird, edge {ec} partially in component{C}")
            #exit(1)
    dists = dict(nx.all_pairs_shortest_path_length(neg_graph))
    #friends - should always be in the same component
    #enemies - should always be in the opposite
    friends = {}
    enemies = {}
    for n in C:
        if n in neg_graph:
            friends[n] = []
            enemies[n] = []
            for other in dists[n]:
                if other != n:
                    #TODO: possibly random starts needed?
                    if dists[n][other] % 2 == 0:
                        friends[n].append(other)
                    else:
                        enemies[n].append(other)
            logger.debug(f"Friends of {n}: {friends[n]}")
            logger.debug(f"Enemies of {n}: {enemies[n]}")

    parts = [set(), set()]
    swap_together = {}

    #TODO: place for randomization
    
    C_nodes = list(C.nodes())
    for n in C_nodes:
        #lets remove non-homologous nodes
        if not n in enemies and n.find("Aux") == -1:
            C.remove_node(n)
            logger.debug(f"Removing node {n}")
        else:
            swap_together[n] = [n]
            start_component = random.randint(0, 1)
            if not n in parts[0] and not n in parts[1]:
                parts[start_component].add(n)
                if n in friends:
                    for f in friends[n]:
                        parts[start_component].add(f)
                    for f in enemies[n]:
                        parts[1 - start_component].add(f)
    
    for n in friends:
        for f in friends[n]:
            swap_together[n].append(f)
        for f in enemies[n]:
            swap_together[n].append(f)    
    
    return parts, swap_together



def optimize_partition(C, parts, swap_together):
    cur_score = score_partition(C, parts)
    #TODO: without randomizations iterations are useless.
    MAX_ITERATIONS = 5    
    not_changed = 0
    all_nodes = []
    for i in range (2):
        for j in parts[i]:
            all_nodes.append(j)
    all_nodes.sort()
    while not_changed < MAX_ITERATIONS:        
        for n in all_nodes:

            score_change = score_swap(C, parts, swap_together[n])
            logger.debug(f"Trying to move {n} and its neighbors {swap_together[n]}, old_score {cur_score}, swap score {score_change}") 
            if score_change > 0:
                upd_score = cur_score + score_change
                for swap_node in swap_together[n]:
                    for j in range (2):
                        if swap_node in parts[j]:
                            parts[j-1].add(swap_node)
                            parts[j].remove(swap_node)
                            break
                not_changed = 0
                logger.debug(f"Moved {n} and its neighbors {swap_together[n]},improving score")
        not_changed += 1

    return parts



#TODO: rewrite, leaving old version for better tests
def create_initial_partition(C, matchGraph):
    p1 = []
    p2 = []
    for n in sorted(C.nodes()):
        if n in matchGraph and len(matchGraph.edges(n)) > 0:
            #TODO: replace with something reasonable! multiple nodes will not work here. Transform each homology component into one node before?.. Anyway, it still should be fixed later with K-L iterations.
            #TODO: distance from starting vertex in match graph % 2
            for ec in matchGraph.edges(n): 
                if matchGraph.edges[ec]['weight'] >= 0:
                    continue
                if ec[0] == n and ec[1] in p1:
                    if n not in p1:
                        p2.append(n)
                elif ec[0] == n and ec[1] in p2:
                    if n not in p2:
                        p1.append(n)
                elif ec[1] == n and ec[0] in p1:
                    if n not in p1:
                        p2.append(n)
                elif ec[1] == n and ec[0] in p2:
                    if n not in p2:
                        p1.append(n)
                elif n not in p1 and n not in p2:
                    if random.random() <= 0.5:
                        p1.append(n)
                    else:
                        p2.append(n)

            if n not in p1 and n not in p2:
                if random.random() <= 0.5:
                    p1.append(n)
                else:
                    p2.append(n)
        else:
            if random.random() <= 0.5:
                p1.append(n)
            else:
                p2.append(n)
#            print("Initial partitions are %s and %s" % (set(p1), set(p2)))
    if len(p1) * len(p2) > 0:
#                    logger.info (f"{len(p1)} {len(p2)} {len(C.nodes())}")
        #Debug for weird crashes
        inter = list(set(p1).intersection(p2))
        missing = set(C.nodes()) - set(p1) - set(p2)
        if len(inter) > 0:
            logger.info (f"Intersecting {inter}")
        if len(missing) > 0:
            logger.info (f"Missing {missing}")
    return [set(p1), set(p2)]

def create_graph_to_phase(current_component, G, matchGraph, hicGraph, uneven_depth, edges, dists):

    C = nx.Graph()
    # rebuild the graph from this component using only edges in the hic graph
    C.add_nodes_from(current_component)

    # first we ignore any nodes that are too short or have too few links.
    short = []

    not_to_big = []
    for n in C.nodes():
        logging.debug(n)
        if G.nodes[n]['coverage'] < MAX_COV:
            not_to_big.append((n, G.nodes[n]))
    local_max_cov = 1.5 * getMedianCov(not_to_big)
    if local_max_cov < 0:
        local_max_cov = MAX_COV
    for n in C.nodes():
        if n not in G:
            sys.stderr.write("Error, got a node not in original graph %s !" % (n))
            sys.exit()
        if G.nodes[n]['coverage'] > local_max_cov and (not uneven_depth):
            logger.debug("While partitoning dropping node %s coverage too high" % (n))
            short.append(n)
        elif not (n in matchGraph) and uneven_depth:
            logger.info("While partitoning dropping node %s uneven coverage and no matches" % (n))
            short.append(n)                
        elif G.nodes[n]['length'] < MIN_LEN:
            collapseOrientedNode(edges, n)
            short.append(n)
        else:
            good = False
            for e in hicGraph.edges(n):
                if (e[0] != n or e[1] != n) and hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT \
                        and (e[0] in C and e[1] in C):
                    good = True
            #Let's phase homologous nodes randomly then
            if n in matchGraph and len(matchGraph.edges(n)) > 0:
                for ec in matchGraph.edges(n): 
                    if matchGraph.edges[ec]['weight'] < 0:
                        good = True
                        break
                            
            if not good:
                logger.info("While partitoning dropping node %s low links count" % (n))
                short.append(n)

    C.remove_nodes_from(short)
    logger.info(f'Currently {C.number_of_nodes()} nodes')
    # Adding some auxilary vertices to allow slightly unbalanced partitions
    # sqrt/2 is quite arbitrary and subject to change
    aux_nodes = int(math.sqrt(C.number_of_nodes() // 2))
    if C.number_of_nodes() <= 3:
        aux_nodes = 0
    for i in range(0, aux_nodes):
        C.add_node("Aux" + str(i))
    for e in hicGraph.edges(current_component):
        # currently only added edges if these nodes are in the component and not matches (homologous) but should allow links to singletons too (to phase disconnected nodes correctly)
        if e[0] in C and e[1] in C and matchGraph.get_edge_data(e[0], e[1]) == None:
            # if edges are too distant in graph, hi-c info is trash
            # using homologous edges when counting distances to help with coverage gaps
            similar_edges = [set(), set()]
            for ind in range(0, 2):
                for match_edge in matchGraph.edges(e[ind]):
                    similar_edges[ind].add(match_edge[0])
                    similar_edges[ind].add(match_edge[1])
                similar_edges[ind].add(e[ind])
#TODO: likely this is not needed anymore since we already added links between homologous edges.
            C.add_edge(e[0], e[1], weight=hicGraph[e[0]][e[1]]['weight'])
            #Possibly we'll need to restore tis for backcomp
#            for e0like in similar_edges[0]:
#                for e1like in similar_edges[1]:
#                    if e0like in dists and e1like in dists[e0like]: #and dists[e0like][e1like] < MAX_GRAPH_DIST + G.nodes[e1like]['length']:
#                        C.add_edge(e[0], e[1], weight=hicGraph[e[0]][e[1]]['weight'])
#                        break
    logger.info(f'Currently {C.number_of_nodes()} in current component')
    logger.debug(C.nodes())
    if C.number_of_nodes() > 1:
        for u, v, w in matchGraph.edges.data("weight"):
            if u in C and v in C:
                if w != None:
                    C.add_edge(u, v, weight=w)
        res = checkXYcomponent(current_component, matchGraph, G, edges)
        if res != [0, 0]:
            sys.stderr.write(f"XY found, {res[0]} {res[1]}, adding fake links\n")
            if res[0] in C and res[1] in C:
                C.add_edge(res[0], res[1], weight=-10 * FIXED_WEIGHT)
                C.add_edge(res[1], res[0], weight=-10 * FIXED_WEIGHT)
        for edge in C.edges:
            logger.debug(f'HIC edge: {edge} {C.edges[edge]}')

    return C

def score_partition(C, partition):
    sum_w = 0
    for part_id  in range(0, 2):
        for i in partition[part_id]:
            for j in partition[part_id]:
                if [i, j] in C.edges():
                    sum_w += C.edges[i, j]['weight']
                    if C.edges[i, j]['weight'] < 0:
                    
                        logger.error(f"Negative edge {i} {j} with weight {C.edges[i, j]['weight']} IN THE SAME PARTITION")
                        exit()
    return sum_w
def score_swap(C, partition, to_swap):
    upd = 0
    for n in to_swap:
        for i in range (2):
            if n in partition[i]:
                for neighbour in C.neighbors(n):
                    if not neighbour in to_swap:
                        if neighbour in partition[i]:
                            upd -= C.edges[n, neighbour]['weight']
                        elif neighbour in partition[1 - i]:
                            upd += C.edges[n, neighbour]['weight']
    return upd

def min_nor_dist(n1, n2, dists):
    res = 10000000000
    for or1 in ['-', '+']:
        for or2 in ['-', '+']:
            or_n1 = n1 + or1
            or_n2 = n2 + or2
            if or_n1 in dists and or_n2 in dists[or_n1]:
                res = min (res, dists[or_n1][or_n2])
    return res

def comnponent_size (current_component, G):
    total_len = 0
    for n in current_component:
        #oriented or not we do not care
        if n in G.nodes():
            total_len += G.nodes[n]['length']
        elif n + '+' in G.nodes():
            total_len += G.nodes[n + '+']['length']
    return total_len

def is_duplication_like(n1, n2, dists, or_G):
    CLOSE_SUSP_HOMOLOGY_DIST = 5000000
    if or_G.nodes[n1 + '+']['length'] > CLOSE_SUSP_HOMOLOGY_DIST and or_G.nodes[n2 + '+']['length'] > CLOSE_SUSP_HOMOLOGY_DIST:
        return
    for or1 in ['-', '+']:
        for or2 in ['-', '+']:
            or_n1 = n1 + or1
            or_n2 = n2 + or2
            if or_n1 in dists and or_n2 in dists[or_n1]:
                dist = (dists[or_n1][or_n2] - or_G.nodes[or_n1]['length'] - or_G.nodes[or_n2]['length']) / 2
                if dist < CLOSE_SUSP_HOMOLOGY_DIST:
                    if or_n2 in dists and or_n1 in dists[or_n2] and dists[or_n2][or_n1] < CLOSE_SUSP_HOMOLOGY_DIST:
                        logger.info(f"Nodes {or_n1} and {or_n2} from haploid component in loop, {dist}")
                    else:
                        logger.info(f"Nodes {or_n1} and {or_n2} one direction close {dist} but not another")
                    return True
    return False

def clear_links(node_pairs, mg):
    for pair in node_pairs:
        mg.matchGraph.remove_edge(pair[0], pair[1])

def remove_pseudo_homology(current_component, or_G, dists, mg):
    #we never remove true homology between really large nodes
    #max 1/3 among all homologies can be cleared 
    MAX_UNDIPLOID_FRACTION = 3

    total_len = comnponent_size(current_component, or_G)
    clear_diploids = []
    total_homology_length = 0
    clear_homology_length = 0
    for n1 in current_component:
        for n2 in current_component:
            if n1 < n2 and mg.isHomologousNodes(n1, n2, True):                
                total_homology_length += mg.getHomologyLength(n1, n2)                                               
                if is_duplication_like(n1, n2, dists, or_G):
                    clear_diploids.append([n1, n2])
                    clear_homology_length += mg.getHomologyLength(n1, n2)
    if (mg.isDiploid(current_component)):
        if clear_homology_length * MAX_UNDIPLOID_FRACTION < total_homology_length:
            if clear_homology_length > 0:
                logger.info(f"Cleaning diploid component {clear_diploids} hom to clear:{clear_homology_length} hom total: {total_homology_length} size total: {total_len}")
            clear_links(clear_diploids, mg)
        else:
            logger.info(f"NOT cleaning diploid component {clear_diploids}, too high fraction. hom to clear:{clear_homology_length} hom total: {total_homology_length} size total: {total_len}")
    else:
        if clear_homology_length > 0:
            logger.info(f"Cleaning haploid component {clear_diploids} hom to clear:{clear_homology_length} hom total: {total_homology_length} size total: {total_len}")
        if clear_homology_length == total_homology_length:
            return "no_diploid"        
    return "diploid"

def process_haploid_component(current_component, G, tsv_file, haploid_component_count):
    haploid_component_count += 1
    haplotype = haploid_component_count % 2
    nodes_to_report = []
    for n in sorted(current_component):
        if G.nodes[n]['length'] > MIN_LEN and G.nodes[n]['coverage'] < MAX_COV:
            nodes_to_report.append(n)
    #TODO: possibly total length condition? Better not assign colors to distal bits...
    if len(nodes_to_report) > 1:
        for n in nodes_to_report:
            report_phased_node(n, haplotype, tsv_file)

def only_aux(C):
    for n in C.nodes():
        if not n.startswith('Aux'):
            return False
    return True

def output_graph_stats(G):
    degrees = [val for (node, val) in G.degree()]
    mean = sum(degrees) / G.number_of_nodes()
    variance = sum([((x - mean) ** 2) for x in degrees]) / G.number_of_nodes()
    res = variance ** 0.5
    logger.info("Loaded a graph with %d nodes and %d edges avg degree %f and stdev %f max is %f" % (
    G.number_of_nodes(), G.number_of_edges(), mean, res, mean + 5 * res))

def run_clustering (graph_gfa, mashmap_sim, hic_byread, output_dir, no_rdna, uneven_depth):
    G = nx.Graph()
    gf.load_indirect_graph(graph_gfa, G)
    output_graph_stats(G)
    delete_set = set()
    largest_component = max(nx.connected_components(G), key=len)
    
    cur_col = 0
    old_colors = {}
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        for e in current_component:
            old_colors[e] = cur_col
        cur_col += 1

#TODO: no_rdna is false, so efficiently we are still removing rDNA component. Should we?
    '''
    if not (no_rdna):
        #Store rDNA component, not to add additional links from matchGraph
        #sys.stderr.write(f"Found an rDNA huge component of {len(largest_component)} edges\n")
        #Here we remove large connected components of short edge, just to exclude rDNA cluster
        delete_set = gf.remove_large_tangles(G, MAX_SHORT_LEN, MAX_SHORT_COMPONENT,MIN_RDNA_COMPONENT, MAX_RDNA_COMPONENT)
    else:
        sys.stderr.write(f"Not using rDNA component removal heuristics\n")
    '''
    filtered_graph = open(os.path.join(output_dir, FILTERED_GFA_FILENAME), 'w')
    tsv_output = os.path.join(output_dir, RESULT_FILENAME)
    tsv_file = open(tsv_output, 'w')
    tsv_file.write("node\tmat\tpat\tmat:pat\tcolor\n")

    translate = open(graph_gfa, 'r')

    for line in translate:
        arr = line.split()
        if arr[0] == "S":
            if not (arr[1] in delete_set):
                filtered_graph.write(line)
        if arr[0] == "L":
            if not(arr[1] in delete_set) and not(arr[3] in delete_set):
                filtered_graph.write(line)
    
    translate.close()

    #loading oriented graph
    #TODO: only place which used it is likely outdated
    #Also we have special function...
    nodelines = []    
    or_G = nx.DiGraph()
    gf.load_direct_graph(graph_gfa, or_G)
    edges = {}
    for l in open(graph_gfa, 'r'):
        parts = l.strip().split('\t')
        if parts[0] == 'S':
            nodelines.append((parts[1], l.strip()))
        elif parts[0] == 'L':
            fromnode = (">" if parts[2] == "+" else "<") + parts[1]
            tonode = (">" if parts[4] == "+" else "<") + parts[3]
   #        edgelines.append((fromnode, tonode, l.strip()))
            if fromnode not in edges:
                edges[fromnode] = set()
            if gf.revnode(tonode) not in edges:
                edges[gf.revnode(tonode)] = set()
            edges[fromnode].add(tonode)
            edges[gf.revnode(tonode)].add(gf.revnode(fromnode))

    #dirty calculation median coverage, considering that most of the phasing is done
    
    med_cov = getMedianCov(G.nodes(data=True))
    MAX_COV = med_cov * 1.5
    
    if (uneven_depth):
        logger.info(f"Will not use coverage based homozygous nodes detection")
    else:
        logger.info(f"Will use coverage based homozygous nodes detection, cutoff: {MAX_COV}")
        
 
    # load hic connections based on mappings, weight corresponds to number of times we see a connection
    hicGraph = gf.loadHiCGraph(hic_byread)
    compressed_file = open(os.path.join(output_dir, HIC_COMPRESSED_FILENAME), 'w')
    max_w = 0
    for node1, node2 in hicGraph.edges():
        if max_w < hicGraph[node1][node2]["weight"]:
            max_w = hicGraph[node1][node2]["weight"]
        compressed_file.write(f'X {node1} {node2} {hicGraph[node1][node2]["weight"]}\n')
    FIXED_WEIGHT = max_w
    logger.info(f'Constant for homologous edges: {-10 * FIXED_WEIGHT}')

    #Adding link between matched edges to include separated sequence to main component

    #TODO: only one of those should be used
    mg = match_graph.MatchGraph(mashmap_sim, G, -10*FIXED_WEIGHT, CLEAR_HOMOLOGY, MIN_ALIGNMENT, logger)
    matchGraph = mg.getMatchGraph()
    
    component_colors = gf.getComponentColors(G)

#reconnecting homologous nodes
    for [v1,v2] in matchGraph.edges():
        if v1 in G.nodes and v2 in G.nodes and matchGraph[v1][v2]['weight'] < 0:    
            if component_colors[v1] != component_colors[v2]:
                logger.info(f"Adding graph link between homologous {v1} {v2}, components {component_colors[v1]} and {component_colors[v2]}")
                G.add_edge(v1, v2)
                


    logger.info(f"Loaded hic info with {hicGraph.number_of_nodes()} nodes and {hicGraph.number_of_edges()} edges")
    compressed_file.close()

    dists = dict(nx.all_pairs_dijkstra_path_length(G, weight=lambda u, v, d: G.nodes[v]['length']))
    logger.info("Distances counted")


    # connected components decomposition and log the IDs and partition each one
    # currently not going to do the right thing on rDNA component
    haploid_component_count = 0
    dists =  dict(nx.all_pairs_dijkstra_path_length(or_G, weight=lambda u, v, d: or_G.edges[u, v]['mid_length']))
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        logger.info("Connected component with %d nodes is: %s" % (len(current_component), current_component))
                
        cleared = remove_pseudo_homology(current_component, or_G, dists, mg)
        #only backward compatibility with KL
        if not mg.isDiploid(current_component):
            haploid_component_count += 1        
        if cleared == "no_diploid":
            #TODO: backwards compatibility not to change distal bits, to_remove
            #if component_size(current_component, G) > 10000000:
            process_haploid_component(current_component, G, tsv_file, haploid_component_count)
            continue

        C = create_graph_to_phase(current_component, G, matchGraph, hicGraph, uneven_depth, edges, dists)

        #TODO: cycle for randomized starts        
        #optimize_partition(C, parts, swap_together)
        
        best_score = - FIXED_WEIGHT * C.number_of_nodes() * C.number_of_nodes()
        best_part = [set(), set()]

        for seed in range(0, KLIN_STARTS):  # iterate on starting partition
            random.seed(seed)
            parts, swap_together = create_initial_partition_noKL(C, matchGraph, seed)
            if only_aux(C):
                logger.info("Only aux nodes left in component, skipping")
                continue            
            part = community.kernighan_lin.kernighan_lin_bisection(C, partition=parts, max_iter=KLIN_ITER, weight='weight', seed=seed)
            logger.debug(f"init_part:\n{sorted(part[0])}\n{sorted(part[1])}")

            sum_w = score_partition(C, part)

            if (sum_w > best_score):
                logger.info(f'Seed {seed} score {sum_w} improved over {best_score}')
                best_part = part
                best_score = sum_w
                logger.debug(f"best_part:\n{sorted(best_part[0])}\n{sorted(best_part[1])}")
        #try to move relatively short edges to fix case of unbalanced haplo sizes (complex repeats on only one haplo)
#        best_part = community.kernighan_lin.kernighan_lin_bisection(C, partition=parts, max_iter=KLIN_ITER, weight='weight', seed=seed)
#        logger.debug(f"parts:\n{sorted(parts[0])}\n{sorted(parts[1])}")
        #logger.debug(f"best_part:\n{sorted(best_part[0])}\n{sorted(best_part[1])}")
        """
        """
        add_part = [set(), set()]
        only_weights = {}
        for n in current_component:
#            if G.nodes[n]['length'] < MIN_LEN and G.nodes[n]['coverage'] < MAX_COV:
            if (not n in C.nodes()) and G.nodes[n]['coverage'] < MAX_COV:

                weights = [0, 0]
                all_weights = [0, 0]
                total_w = 0
                for e in hicGraph.edges(n):
                    if (e[0] != n or e[1] != n) and (e[0] in current_component and e[1] in current_component):
                        #logger.debug(f'edge {e} to short {hicGraph[e[0]][e[1]]["weight"]}')
                        for ind in range(0, 2):
                            if e[0] in best_part[ind] or e[1] in best_part[ind]:
                                if hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT:
                                    logger.debug(f'{e} edge to partition {ind}!')
                                    weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                                all_weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                        total_w += hicGraph[e[0]][e[1]]["weight"]
                if weights[0] > SIGNIFICANT_MAJORITY * weights[1]:
                    add_part[0].add(n)
                    logger.info(f"Edge {n} assigned color 0")
                elif weights[1] > SIGNIFICANT_MAJORITY * weights[0]:
                    add_part[1].add(n)
                    logger.debug(f"Edge {n} assigned color 1")
                else:
                    logger.debug(f"Edge {n} weights {weights[0]} and {weights[1]} coverage {G.nodes[n]['coverage']}")
                    logger.debug(f"Edge {n} real weights {all_weights[0]} and {all_weights[1]} coverage {G.nodes[n]['coverage']} total weight {total_w}")
                    if (all_weights[0] > SIGNIFICANT_MAJORITY * all_weights[1] or all_weights[1] > SIGNIFICANT_MAJORITY * all_weights[0]) and (all_weights[0] > MIN_WEIGHT or all_weights[1]> MIN_WEIGHT):
                        logger.debug(f"Edge {n} real weights allows to make a decision forbidden by weights!")
                        if G.nodes[n]['length'] < MIN_LEN:
                            only_weights[n] = all_weights
                        else:
                            if all_weights[0] > all_weights[1]:
                                add_part[0].add(n)
                                logger.debug(f"Edge {n} assigned color 0 by additional all_weights")        
                            else:
                                add_part[1].add(n)
                                logger.debug(f"Edge {n} assigned color 1 by additional all_weights")
            elif G.nodes[n]['length'] < MIN_LEN:
                logger.debug(f"Edge {n} did not pass coverage check")

        logger.info(f'Added {len(add_part[0]) + len (add_part[1])} short edges to bipartition')
        logger.info(f'RES\t{add_part}')
        best_part[0].update(add_part[0])
        best_part[1].update(add_part[1])
        for ind in range (0, 2):
            for contig in sorted(best_part[ind]):
                #filtering out Aux nodes
                if contig in current_component:
                    report_phased_node(contig, ind, tsv_file)
        for contig in sorted(only_weights.keys()):
            tsv_file.write(f'{contig}\t{only_weights[contig][0]}\t{only_weights[contig][1]}\t{only_weights[contig][0]}:{only_weights[contig][1]}\t#88FF88\n')


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(f'Usage: {sys.argv[0]} graph.gfa homologous_nodes.matches hic_byread output_dir, no_rdna, uneven_depth')
        exit()
    run_clustering(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
