import networkx as nx
import sys


def get_component_length(G, component):
    res = 0
    for node in component:
        res += G.nodes[node]['length']
    return res

#Functions for indirect graph processing, nodes - "utig1-234"

#TODO: here should be a better check, with either rDNA nodes or diploid structure. After scaff refactor
def remove_large_tangles(G, MAX_LEN, MAX_SHORT_COMPONENT, MIN_RDNA_COMPONENT, MAX_RDNA_COMPONENT):
    shorts = set()
    for node in G.nodes():
        if G.nodes[node]['length'] < MAX_LEN:
            shorts.add(node)
    sh_G = G.subgraph(shorts)
    nodes_deleted = 0
    components_deleted = 0
    to_delete = []
    #largest tangle will be always removed (unless --no-rdna-tangle used which completely disable tangle removal heuristics)
    First = True

    for comp in sorted(nx.connected_components(sh_G), key=len, reverse=True):
        if len(comp) > MAX_SHORT_COMPONENT:
#lets check whether component looks similar to rdna 
            first_edge = list(comp)[0]
            sys.stderr.write(f'Checking component with first edge {first_edge}\n')
            first_edge_big_comp = nx.node_connected_component(G, first_edge)
            not_deleting = set( )
            for e in first_edge_big_comp:
                if not e in comp:
                    not_deleting.add(e)
            subgraph_comp = G.subgraph(not_deleting)
            new_connected = nx.connected_components(subgraph_comp)
            good_new_connected = False
            for new_comp in new_connected:
                tlen = get_component_length(G, new_comp)
                if tlen < MAX_RDNA_COMPONENT and tlen > MIN_RDNA_COMPONENT:
                    good_new_connected = True                     
                    sys.stderr.write(f'Component  size {tlen}, small enough\n')               
                    break
            if First or good_new_connected:
                First = False                   
                components_deleted += 1
                for e in comp:
                    nodes_deleted += 1
                    to_delete.append(e)
            else:
                sys.stderr.write(f'Not deleting component {comp}, nothing similar to rDNA\n')
    G.remove_nodes_from(to_delete)
    sys.stderr.write(f'Removed {components_deleted} short nodes components and {nodes_deleted} short nodes. New '
                     f'number of nodes {G.number_of_nodes()}\n')
    return set(to_delete)

def nodes_in_tangles(G, MAX_LEN, MIN_TANGLE_SIZE):
    shorts = set()
    for node in G.nodes():
        if G.nodes[node]['length'] < MAX_LEN:
            shorts.add(node)
    res = set()
    sh_G = G.subgraph(shorts)
    for comp in nx.connected_components(sh_G):
        if len(comp) > MIN_TANGLE_SIZE:
            for n in comp:
                res.add(n)
    return res

def tsv2gaf(tsv_path):
    res = ""
    arr = tsv_path.strip().split(',')
    for node in arr:
        if node.endswith("+"):
            res += ">" +node[:-1]
        elif node.endswith("-"):
            res += "<" +node[:-1]
        else:
            res += node
    return res
        
#TODO: get rid of code dup
def load_indirect_graph(gfa_file, G):

    for line in open(gfa_file, 'r'):
        if "#" in line:
            continue
        line = line.strip().split()

        if line[0] == "S":
            #noseq graphs
            ls = len(line[2])
            cov = 0
            for i in range(3, len(line)):
                spl_tag = line[i].split(":")
                if spl_tag[0] == "LN":
                    ls = int(spl_tag[2])
                if spl_tag[0] == "ll":
                    cov = float(spl_tag[2])                           
            G.add_node(line[1], length=int(ls), coverage=float(cov))
            
    #L and S lines can be mixed
    for line in open(gfa_file, 'r'):
        if "#" in line:
            continue
        line = line.strip().split()       
        if line[0] == "L":
            if line[1] not in G or line[3] not in G:
                sys.stderr.write("Error while graph loading; link between nodes not in graph:%s" % (line))
                sys.exit(1)
            overlap = 0
            if len(line) > 5 and isinstance(line[5][:-1], int):
                overlap = int(line[5][:-1])
            #Double distance between centers of nodes
            G.add_edge(line[1], line[3], mid_length = G.nodes[line[1]]['length'] + G.nodes[line[3]]['length'] - 2*overlap)

def load_direct_graph(gfa_file, G):
    rc = {"+": "-", "-": "+"}
    
    
    for line in open(gfa_file, 'r'):
        if "#" in line:
            continue
        line = line.strip().split()

        if line[0] == "S":
            #noseq graphs
            ls = len(line[2])
            cov = 0
            for i in range(3, len(line)):
                spl_tag = line[i].split(":")
                if spl_tag[0] == "LN":
                    ls = int(spl_tag[2])
                if spl_tag[0] == "ll":
                    cov = float(spl_tag[2])
            for suff in ["+", "-"]:
                G.add_node(line[1]+suff, length=int(ls), coverage=float(cov))                        
            
    #L and S lines can be mixed
    for line in open(gfa_file, 'r'):
        line = line.strip().split()       
        if line[0] == "L":
            if line[1]+"+" not in G or line[3]+"+" not in G:
                sys.stderr.write("Error while graph loading; link between nodes not in graph:%s" % (line))
                sys.exit(1)
            overlap = 0
            if len(line) > 5 and isinstance(line[5][:-1], int):
                overlap = int(line[5][:-1])
            #Double distance between centers of nodes
            G.add_edge(line[1] + line[2], line[3] + line[4], mid_length = G.nodes[line[1]+line[2]]['length'] + G.nodes[line[3]+line[4]]['length'] - 2*overlap)
            G.add_edge(line[3] + rc[line[4]], line[1] + rc[line[2]], mid_length = G.nodes[line[1]+line[2]]['length'] + G.nodes[line[3]+line[4]]['length'] - 2*overlap)

def loadHiCGraph(hic_byread_file):
    hicGraph = nx.Graph()
    hic_file = open(hic_byread_file, 'r')
    for line in hic_file:
        if "#" in line:
            continue
        line = line.strip().split()

        if len(line) < 3:
            continue

        if line[1] == line[2]:
            continue

        hicGraph.add_node(line[1])
        hicGraph.add_node(line[2])
        if len(line) > 3:
            add_w = int(line[3])
        else:
            add_w = 1
        w = hicGraph.get_edge_data(line[1], line[2], 0)
        if w == 0:
            hicGraph.add_edge(line[1], line[2], weight=add_w)
        else:
            w = w['weight'] + add_w
            hicGraph[line[1]][line[2]]['weight'] = w
    return hicGraph

def getComponentColors(G):
    component_colors = {}
    current_color = 0

    for current_component in sorted(nx.connected_components(G), key=len, reverse=True):
        for e in current_component:
            component_colors[e] = current_color
        current_color += 1
    return component_colors

# homology_weight - large negative something, homology_len - minimal length of homology to be taken in account. 
# Some shorter ones can still sometimes be used if they are in regular bulge_like structure
def loadMatchGraph(mashmap_sim, G, homology_weight, homology_len):
    # from node to [best_match, weight, second_best_weight], second best to other node!
    # Possibly will have to unite different matches between same node pairs splitted by mashmap
    best_match = {}
    mashmap_weights = {}
    matchGraph = nx.Graph()
    component_colors = getComponentColors(G)

    for line in open(mashmap_sim, 'r'):
        if "#" in line:
            continue
        line = line.strip().split()
        if len(line) < 3:
            continue
        if line[0] == line[1]:
            continue
        if not line[0] in best_match:
            best_match[line[0]] = [line[1], int(line[2]), 1]
        if not line[1] in best_match:
            best_match[line[1]] = [line[0], int(line[2]), 1]
        if int(line[2]) > best_match[line[0]][1]:
            if best_match[line[0]][0] == line[1]:
                # not updating second best match
                best_match[line[0]] = [line[1], int(line[2]), best_match[line[0]][2]]
            else:
                best_match[line[0]] = [line[1], int(line[2]), best_match[line[0]][1]]

        if int(line[2]) > best_match[line[1]][1]:
            if best_match[line[1]][0] == line[0]:
                best_match[line[1]] = [line[0], int(line[2]), best_match[line[1]][2]]
            else:
                best_match[line[1]] = [line[0], int(line[2]), best_match[line[1]][1]]

        # only debug purpose
        upd_weight = int(line[2])
        if line[0] + line[1] in mashmap_weights:
            upd_weight = max(upd_weight, mashmap_weights[line[0] + line[1]])
        mashmap_weights[line[0] + line[1]] = upd_weight
        mashmap_weights[line[1] + line[0]] = upd_weight

        if int(line[2]) > homology_len:
            color_sum = 0
            for id in range(0, 2):
                if line[id] in component_colors:
                    color_sum += 1
            if color_sum != 2:
                sys.stderr.write(f"Weird deleted node pair {line[0]} {line[1]}\n")
                continue

            matchGraph.add_edge(line[0], line[1])

    # Now let's add links for shorter nodes
    for line in open(mashmap_sim, 'r'):
        if "#" in line:
            continue
        line = line.strip().split()
        if len(line) < 3:
            continue
        neighbours = [set(), set()]
        for id in range(0, 2):
            for edge in G.edges(line[id]):
                neighbours[id].add(edge[0])
                neighbours[id].add(edge[1])
        common = neighbours[0].intersection(neighbours[1])
        if len(common) > 0:
            # we are in somethiing bulge-like
            len_cutoff = min(G.nodes[line[0]]['length'], G.nodes[line[1]]['length'], homology_len) / 2
        # possibly check whether we have something already phased nearby?
            if int(line[2]) > len_cutoff:
                matchGraph.add_edge(line[0], line[1])



    for ec in matchGraph.edges():
        # while we build the initial partition give a big bonus edge for putting the homologous nodes into different partitions             
        # Adding an edge that already exists updates the edge data (in networkX graph)
        # At least one in the pair is the best similarity match for other (and second best is not too close to the best)
        #,consecutive edges never assigned to be homologous
        if (not G.has_edge(ec[0], ec[1])) and ((ec[0] in best_match and ec[1] == best_match[ec[0]][0] and best_match[ec[0]][2]/best_match[ec[0]][1]< 0.8) or (ec[1] in best_match and ec[0] == best_match[ec[1]][0] and best_match[ec[1]][2]/best_match[ec[1]][1]< 0.8)):
            matchGraph.add_edge(ec[0], ec[1], weight = homology_weight)
        else:
        #not really look like homologous node pair but still suspicious, lets just wipe the hi-c links but not prioritize splitting them to different partitions.                                
            matchGraph.add_edge(ec[0], ec[1], weight=0)
            
    #Just debug information about homology link structure
    debug_neighbors = {}              
    for debug_node in matchGraph.nodes():
        debug_neighbors[debug_node] = set()
    for debug_node in matchGraph.nodes():
        for edge_pair in matchGraph.edges(debug_node):  
            if matchGraph.edges[edge_pair]['weight'] < 0:                      
                debug_neighbors[edge_pair[0]].add(edge_pair[1])
                debug_neighbors[edge_pair[1]].add(edge_pair[0]) 
                if edge_pair[0] == edge_pair[1]:
                    sys.stderr.write(f"loop {edge_pair} {matchGraph.edges[edge_pair]}\n")
    for node1 in debug_neighbors:
        for node2 in debug_neighbors[node1]:
            for node3 in debug_neighbors[node2]:
                if node3 in debug_neighbors[node1] and node1 < node2 and node2 < node3:
                    sys.stderr.write(f"Found a triangle {node1} {node2} {node3}\n")
                if node1 < node3:
                    sys.stderr.write (f"Multiple homology from {node2}, {node1} with {mashmap_weights[node1+node2]} and {node3} with {mashmap_weights[node3+node2]}\n")
    sys.stderr.write("Loaded match info with %d nodes and %d edges\n" % (matchGraph.number_of_nodes(), matchGraph.number_of_edges()))
    return matchGraph

