import sys
import random
import networkx as nx
import math
import os
import graph_functions
import copy
from networkx.algorithms import community
import graph_functions as gf
from scaffolding import match_graph, logger_wrap

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
        sys.stderr.write("MULTIPLE CANDIDATES FOR FAKE HOMOLOGY IN XY\n")
        for pair in fake_hom:
            sys.stderr.write(f"{pair[0]} --- {pair[1]}\n")
        return [0, 0]
    elif len(fake_hom) == 1:
        sys.stderr.write(f"Adding FAKE HOMOLOGY between {fake_hom[0][0]} --- {fake_hom[0][1]}\n")   
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
        print (f'Node {node[0]} coverage {node[1]["coverage"]} length {node[1]["length"]}')
        sum_l += node[1]['length']
        if 2*sum_l > total_length:
            med_cov = node[1]['coverage']
            print (f'Median coverage is {med_cov}\n')
            break
    return med_cov

def run_clustering (graph_gfa, mashmap_sim, hic_byread, output_dir, no_rdna, uneven_depth):
    #TODO: move constants to some more relevant place
    MIN_LEN = 200000  # best result so far with 200000

    #likely do not need
    #MAX_GRAPH_DIST = 100000000  # hi-c links over 10M are believed to be useless

    KLIN_STARTS = 1000  # number of different starts of kernighan lin
    KLIN_ITER = 10000  # number of iterations inside kernighan lin

    MAX_SHORT_COMPONENT = 100 # we remove from consideration each connected compoents of edges < MIN_SHORT_LEN that is larger than MAX_SHORT_COMPONENT
    MAX_SHORT_LEN=200000 #we remove from consideration each connected compoents of edges < MIN_SHORT_LEN that is larger than
    MIN_WEIGHT = 10 # Ignore edges with few links
    SIGNIFICANT_MAJORITY = 2.5

    RESULT_FILENAME = "hicverkko.colors.tsv"
    FILTERED_GFA_FILENAME = "filtered.gfa"
    HIC_COMPRESSED_FILENAME = "hic.byread.compressed"
    LOGGING_FILENAME = "hicverkko.log"

    CLEAR_HOMOLOGY = 500000 #for clear bulge like structures this limmit is decreased
    MIN_ALIGNMENT = 100000 #smaller alingments will be filtered out
    MAX_COV = 100  # tempora# ry coverage cutoff, currently replaced by median coverage from gfa
    FIXED_WEIGHT = 100000  # best result so far with 100000 #currently replaced with max pairwise weight among datasets

    MAX_RDNA_COMPONENT = 10000000 # maximal size of rDNA component, used for filtering out rDNA cluster only
    MIN_RDNA_COMPONENT = 500000 
    logger = logger_wrap.initLogger("phasing.log")

    # load the assembly gfa
    G = nx.Graph()
    logging_f = open (os.path.join(output_dir, LOGGING_FILENAME), 'w')
    graph_functions.load_indirect_graph(graph_gfa, G)

    degrees = [val for (node, val) in G.degree()]
    mean = sum(degrees) / G.number_of_nodes()
    variance = sum([((x - mean) ** 2) for x in degrees]) / G.number_of_nodes()
    res = variance ** 0.5
    sys.stderr.write("Loaded a graph with %d nodes and %d edges avg degree %f and stdev %f max is %f\n" % (
    G.number_of_nodes(), G.number_of_edges(), mean, res, mean + 5 * res))

    delete_set = set()
    largest_component = max(nx.connected_components(G), key=len)
    cur_col = 0
    old_colors = {}
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        for e in current_component:
            old_colors[e] = cur_col
        cur_col += 1

#TODO: no_rdna is false, so efficiently we are still removing rDNA component. Should we?
    if not (no_rdna):
        #Store rDNA component, not to add additional links from matchGraph
        #sys.stderr.write(f"Found an rDNA huge component of {len(largest_component)} edges\n")
        #Here we remove large connected components of short edge, just to exclude rDNA cluster
        delete_set = graph_functions.remove_large_tangles(G, MAX_SHORT_LEN, MAX_SHORT_COMPONENT,MIN_RDNA_COMPONENT, MAX_RDNA_COMPONENT)
    else:
        sys.stderr.write(f"Not using rDNA component removal heuristics\n")
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


#let's find that nodes that have multiple extension
    multiple_ext = set()
    for node in edges:
        if len(node) > 1:
            multiple_ext.add(node)
    #dirty calculation median coverage, considering that most of the phasing is done
    
    sorted_nodes = sorted(G.nodes(data=True), key=lambda node: node[1]['length'])
    med_cov = getMedianCov(G.nodes(data=True))
    MAX_COV = med_cov * 1.5
    
    if (uneven_depth):
        sys.stderr.write(f"Will not use coverage based homozygous nodes detection\n")
    else:
        sys.stderr.write(f"Will use coverage based homozygous nodes detection, cutoff: {MAX_COV}\n")
        
 
    # load hic connections based on mappings, weight corresponds to number of times we see a connection
    hicGraph = graph_functions.loadHiCGraph(hic_byread)
    compressed_file = open(os.path.join(output_dir, HIC_COMPRESSED_FILENAME), 'w')
    max_w = 0
    for node1, node2 in hicGraph.edges():
        if max_w < hicGraph[node1][node2]["weight"]:
            max_w = hicGraph[node1][node2]["weight"]
        compressed_file.write(f'X {node1} {node2} {hicGraph[node1][node2]["weight"]}\n')
    FIXED_WEIGHT = max_w
    sys.stderr.write(f'Constant for neighbouring edges set to be  {FIXED_WEIGHT} (but not used), for homologous edges {-10 * FIXED_WEIGHT} \n')

    #Adding link between matched edges to include separated sequence to main component

    #TODO: only one of those should be used
    mg = match_graph.MatchGraph(mashmap_sim, G, -10*FIXED_WEIGHT, CLEAR_HOMOLOGY, MIN_ALIGNMENT, logger)
    matchGraph = mg.getMatchGraph()
    
    component_colors = graph_functions.getComponentColors(G)

#reconnecting homologous nodes
    for [v1,v2] in matchGraph.edges():
        if v1 in G.nodes and v2 in G.nodes and matchGraph[v1][v2]['weight'] < 0:    
            if component_colors[v1] != component_colors[v2]:
                logging_f.write(f"Adding graph link between homologous {v1} {v2}, components {component_colors[v1]} and {component_colors[v2]}\n")
                G.add_edge(v1, v2)
                


    sys.stderr.write("Loaded hic info with %d nodes and %d edges\n" % (hicGraph.number_of_nodes(), hicGraph.number_of_edges()))
    compressed_file.close()

    dists = dict(nx.all_pairs_dijkstra_path_length(G, weight=lambda u, v, d: G.nodes[v]['length']))
    sys.stderr.write("Distances counted\n")


    # connected components decomposition and log the IDs and partition each one
    # currently not going to do the right thing on rDNA component
    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        logging_f.write("Connected component with %d nodes is: %s\n" % (len(current_component), current_component))

        total_l = 0
        for n in current_component:
            total_l += G.nodes[n]['length']
        #TODO: currently just a debug message!
        # we likely need to skip phasing for non-diploids and assign the same color for the component, but this can cause troubles for short arms until we stop removing rDNA tangle.
        if total_l > 1000000 and not mg.isDiploid(current_component):
            logging_f.write(f"Component is not diploid!\n")
            

        C = nx.Graph()
        # rebuild the graph from this component using only edges in the hic graph
        C.add_nodes_from(current_component)

        # first we ignore any nodes that are too short or have too few links.
        short = []

        not_to_big = []
        for n in C.nodes():
            print (n)
            if G.nodes[n]['coverage'] < MAX_COV:
                not_to_big.append((n, G.nodes[n]))
        local_max_cov = 1.5 * getMedianCov(not_to_big)
        if local_max_cov < 0:
            local_max_cov = MAX_COV
        for n in C.nodes():
            if n not in G:
                sys.stderr.write("Error got a node not in original graph %s !" % (n))
                sys.exit()
            if G.nodes[n]['coverage'] > local_max_cov and (not uneven_depth):
                logging_f.write("While partitoning dropping node %s coverage too high\n" % (n))
                short.append(n)
            elif not (n in matchGraph) and uneven_depth:
                logging_f.write("While partitoning dropping node %s uneven coverage and no matches\n" % (n))
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
                    logging_f.write("While partitoning dropping node %s low links count\n" % (n))
                    short.append(n)

        C.remove_nodes_from(short)
        logging_f.write(f'Currently {C.number_of_nodes()} nodes\n')
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
                for e0like in similar_edges[0]:
                    for e1like in similar_edges[1]:
                        if e0like in dists and e1like in dists[e0like]: #and dists[e0like][e1like] < MAX_GRAPH_DIST + G.nodes[e1like]['length']:
                            C.add_edge(e[0], e[1], weight=hicGraph[e[0]][e[1]]['weight'])
                            break
        logging_f.write(f'Currently {C.number_of_nodes()} in current component\n')

        
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
                logging_f.write(f'HIC edge: {edge} {C.edges[edge]}\n')


            best_score = FIXED_WEIGHT * C.number_of_nodes() * C.number_of_nodes()




            for seed in range(0, KLIN_STARTS):  # iterate on starting partition
                random.seed(seed)
                p1 = []
                p2 = []
                for n in sorted(C.nodes()):
                    if n in matchGraph and len(matchGraph.edges(n)) > 0:
                        #TODO: replace with something reasonable! multiple nodes will not work here. Transform each homology component into one node before?.. Anyway, it still should be fixed later with K-L iterations.
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
#                    logging_f.write (f"{len(p1)} {len(p2)} {len(C.nodes())}\n")
                    #Debug for weird crashes
                    inter = list(set(p1).intersection(p2))
                    missing = set(C.nodes()) - set(p1) - set(p2)
                    if len(inter) > 0:
                        logging_f.write (f"Intersecting {inter}\n")
                    if len(missing) > 0:
                        logging_f.write (f"Missing {missing}\n")

                    part = community.kernighan_lin.kernighan_lin_bisection(C, partition=[set(p1), set(p2)], max_iter=KLIN_ITER,
                                                                           weight='weight', seed=seed)
                    sum_w = 0
                    for i in part[0]:
                        for j in part[1]:
                            if [i, j] in C.edges():
                                # print (f'{i} {j} edge   {C.edges[i,j]}')
                                sum_w += C.edges[i, j]['weight']

                    if (sum_w < best_score):
    #lets forbid Aux only components
                        if check_non_empty(part[0], G) and check_non_empty(part[1], G):
                            logging_f.write(f'Seed {seed} score {sum_w} improved over {best_score}\n')
                            best_part = part
                            best_score = sum_w
            #try to move relatively short edges to fix case of unbalanced haplo sizes (complex repeats on only one haplo)
            fixUnbalanced(best_part, C, G)
            logging_f.write(f'RES\t{best_part}\n')

            add_part = [set(), set()]
            only_weights = {}
            for n in current_component:
                if G.nodes[n]['length'] < MIN_LEN and G.nodes[n]['coverage'] < MAX_COV:
                    weights = [0, 0]
                    all_weights = [0, 0]
                    total_w = 0
                    for e in hicGraph.edges(n):
                        if (e[0] != n or e[1] != n) and (e[0] in current_component and e[1] in current_component):
                            logging_f.write(f'DEBUG: edge {e} to short {hicGraph[e[0]][e[1]]["weight"]}\n')
                            for ind in range(0, 2):
                                if e[0] in best_part[ind] or e[1] in best_part[ind]:
                                    if hicGraph[e[0]][e[1]]["weight"] > MIN_WEIGHT:
                                        logging_f.write(f'DEBUG: {e} edge to partition {ind}!\n')
                                        weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                                    all_weights[ind] += hicGraph[e[0]][e[1]]["weight"]
                            total_w += hicGraph[e[0]][e[1]]["weight"]
                    if weights[0] > SIGNIFICANT_MAJORITY * weights[1]:
                        add_part[0].add(n)
                        logging_f.write(f"Edge {n} assigned color 0\n")
                    elif weights[1] > SIGNIFICANT_MAJORITY * weights[0]:
                        add_part[1].add(n)
                        logging_f.write(f"Edge {n} assigned color 1\n")
                    else:
                        logging_f.write(f"Edge {n} weights {weights[0]} and {weights[1]} coverage {G.nodes[n]['coverage']}\n")
                        logging_f.write(f"Edge {n} real weights {all_weights[0]} and {all_weights[1]} coverage {G.nodes[n]['coverage']} total weight {total_w}\n")
                        if (all_weights[0] > SIGNIFICANT_MAJORITY * all_weights[1] or all_weights[1] > SIGNIFICANT_MAJORITY * all_weights[0]) and (all_weights[0] > MIN_WEIGHT or all_weights[1]> MIN_WEIGHT):
                            logging_f.write(f"Edge {n} real weights allows to make a decision forbidden by weights!\n")
                            only_weights[n] = all_weights
                elif G.nodes[n]['length'] < MIN_LEN:
                    logging_f.write(f"Edge {n} did not pass coverage check\n")

            logging_f.write(f'Added {len(add_part[0]) + len (add_part[1])} short edges to bipartition\n')
            logging_f.write(f'RES\t{add_part}\n')
            best_part[0].update(add_part[0])
            best_part[1].update(add_part[1])
            right = False
            for ind in range (0, 2):
                for contig in sorted(best_part[ind]):
                    if contig in current_component:
                        if ind == 1:
                            tsv_file.write(f'{contig}\t0\t100000\t0:100000\t#8888FF\n')
                        else:
                            tsv_file.write(f'{contig}\t100000\t0\t100000:0\t#FF8888\n')
            for contig in sorted(only_weights.keys()):
                tsv_file.write(f'{contig}\t{only_weights[contig][0]}\t{only_weights[contig][1]}\t{only_weights[contig][0]}:{only_weights[contig][1]}\t#88FF88\n')


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(f'Usage: {sys.argv[0]} graph.gfa homologous_nodes.matches hic_byread output_dir, no_rdna, uneven_depth')
        exit()
    run_clustering(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
