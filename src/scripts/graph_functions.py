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

# https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
trans = str.maketrans("ACTG", "TGAC")
def revcomp(s):
    return s.translate(trans)[::-1]

def str2bool(v):
    return v.lower() in ("yes", "true", "t", "1")

def revnode(n):
    assert len(n) >= 2
    assert n[0] == "<" or n[0] == ">"
    return (">" if n[0] == "<" else "<") + n[1:]

def find(parent, key):
    assert key in parent
    while parent[key] != parent[parent[key]]:
        parent[key] = parent[parent[key]]
    return parent[key]

def merge(parent, left, right):
    left = find(parent, left)
    right = find(parent, right)
    parent[right] = left

def canon(left, right):
    fwstr = left + right
    bwstr = revnode(right) + revnode(left)
    if bwstr < fwstr: return (revnode(right), revnode(left))
    return (left, right)

def canontip(left, right):
    fwstr = left + right
    bwstr = right + left
    if bwstr < fwstr: return (right, left)
    return (left, right)

def getone(s):
    assert len(s) >= 1
    for i in s:
        return i

def pathstr(p):
    return "".join(p)

def iterate_deterministic(l, end = ""):
    tmp = list(l)
    tmp.sort()

    if end != "" and end in tmp:
        tmp.remove(end)
        tmp.insert(0, end)
    for x in tmp:
        yield x

# https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
# iterative so we don't hit max recursion depth
def strong_connect_iterative(start, i, index, lowlink, on_stack, result, edges):
    stack = []
    stack.append((0, start, 0))
    S = []
    while len(stack) > 0:
        (state, node, neighbor_i) = stack[-1]
        stack.pop()
        if state == 0:
            assert node not in on_stack
            assert node not in index
            assert node not in lowlink
            index[node] = i
            lowlink[node] = i
            i += 1
            S.append(node)
            on_stack.add(node)
            stack.append((1, node, 0))
        elif state == 1:
            if node in edges and neighbor_i < len(edges[node]):
                neighbor = edges[node][neighbor_i][0]
                if neighbor not in index:
                    stack.append((2, node, neighbor_i))
                    stack.append((0, neighbor, 0))
                    continue
                elif neighbor in on_stack:
                    assert neighbor in index
                    assert node in lowlink
                    lowlink[node] = min(lowlink[node], index[neighbor])
                neighbor_i += 1
            if node in edges and neighbor_i < len(edges[node]):
                stack.append((1, node, neighbor_i))
            else:
                stack.append((3, node, 0))
        elif state == 2:
            neighbor = edges[node][neighbor_i][0]
            assert neighbor in lowlink
            lowlink[node] = min(lowlink[node], lowlink[neighbor])
            neighbor_i += 1
            stack.append((1, node, neighbor_i))
        elif state == 3:
            assert node in lowlink
            assert node in index
            if lowlink[node] == index[node]:
                result.append([])
                stacknode = ""
                while stacknode != node:
                    assert len(S) > 0
                    stacknode = S[-1]
                    S.pop()
                    assert stacknode in on_stack
                    on_stack.remove(stacknode)
                    result[-1].append(stacknode)
        else:
            assert False
    assert len(S) == 0
    return i

# https://en.wikipedia.org/wiki/Tarjan%27s_strongly_connected_components_algorithm
def topological_sort(nodelens, edges):
    index = {}
    lowlink = {}
    on_stack = set()
    result = []
    i = 0
    for node in nodelens:
        if ">" + node not in index:
            old_i = i
            i = strong_connect_iterative(">" + node, i, index, lowlink, on_stack, result, edges)
            assert i > old_i
        if "<" + node not in index:
            old_i = i
            i = strong_connect_iterative("<" + node, i, index, lowlink, on_stack, result, edges)
            assert i > old_i
    result = result[::-1]
    belongs_to_component = {}
    for i in range(0, len(result)):
        assert len(result[i]) >= 1
        for node in result[i]:
            assert node not in belongs_to_component
            belongs_to_component[node] = i
    for node in nodelens:
        assert ">" + node in belongs_to_component
        assert "<" + node in belongs_to_component
    for node in belongs_to_component:
        if node in edges:
            for edge in edges[node]:
                assert belongs_to_component[edge[0]] >= belongs_to_component[node]
    return (result, belongs_to_component)

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


#classes for matchGraph construction
class HomologyInfo:
    def __init__(self, node1, node2, len1, len2):
        self.nodes = [node1, node2]
        self.len = [len1, len2]
        self.covered = [0, 0]
        #two lists for each nodes, each list constains [st_i, 1] and [end_i, -11]
        self.intervals = [[],[]]

    def addInterval(self, intervals):
        for i in range(0, 2):
            self.intervals[i].append([intervals[i][0],1])
            self.intervals[i].append([intervals[i][1],-1])

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
    
    def getMinCovered(self):
#in weird case matches can be larger than seq, avoiding

        return min(self.covered[0], self.covered[1], self.len[0], self.len[1])
    
    def getMinLength(self):
        return min(self.len[0], self.len[1])

class HomologyStorage:
    #{node1: {node2: HomologyInfo(node_1, node_2)}}
    def __init__(self):
        self.homologies = {}

    def addHomology(self, node1, node2, len1, len2, intervals):
        if not node1 in self.homologies:
            self.homologies[node1] = {}
        if not node2 in self.homologies[node1]:
            self.homologies[node1][node2] = HomologyInfo(node1, node2, len1, len2)
        self.homologies[node1][node2].addInterval(intervals)
    
    def fillCoverage(self):
        for node1 in self.homologies:
            for node2 in self.homologies[node1]:    
                self.homologies[node1][node2].fillCoverage()


# homology_weight - large negative something, min_big_homology - minimal length of homology to be taken in account. 
# Some shorter ones can still sometimes be used if they are in regular bulge_like structure
def loadMatchGraph(mashmap_sim, G, homology_weight, min_big_homology, min_alignment):
    matchGraph = nx.Graph()
    hom_storage = HomologyStorage()
    for line in open(mashmap_sim, 'r'):
        arr = line.strip().split()
        if len(arr) < 11:
            continue
        if int(arr[10]) < min_alignment:
            continue
        #self mapping
        if arr[0] == arr[5]:
            continue
        if len(arr[0]) < 3:
            print (line)
            exit()
        #utig4-0 2145330 0       990000  +       utig4-0 2145330 12      994065  37      994053  51      id:f:0.999992   kc:f:0.874893
        hom_storage.addHomology(arr[0], arr[5], int(arr[1]), int(arr[6]), [[int(arr[2]), int(arr[3])], [int(arr[7]), int(arr[8])]])
    hom_storage.fillCoverage()

    for node1 in hom_storage.homologies:
        for node2 in hom_storage.homologies[node1]:
            #we deleted some nodes after mashmap
            #sys.stderr.write(f"Checking {node1} {node2} {hom_storage.homologies[node1][node2].getMinCovered()} {min_big_homology}\n")
            if node1 in G.nodes() and node2 in G.nodes():
                cur_homology = hom_storage.homologies[node1][node2].getMinCovered()
                if cur_homology > min_big_homology:
                    matchGraph.add_edge(node1, node2, homology_len = cur_homology)
                else:
                    #less strict condition for bulge-like structure
                    nodes = [node1, node2]
                    neighbours = [set(), set()]
                    for i in range(0, 2):
                        for adj_node in G.neighbors(nodes[i]):
                            neighbours[i].add(adj_node)
                    #common = neighbours[0].intersection(neighbours[1])
                    len_cutoff = min(hom_storage.homologies[node1][node2].getMinLength(), min_big_homology)/2
                    # possibly check whether we have something already phased nearby?
                    # very strict bulge-like condition
                    if cur_homology > len_cutoff and neighbours[0] == neighbours[1] and len(neighbours[0]) > 0:
                        sys.stderr.write(f"Adding bulge-like {node1} {node2} {cur_homology} {len_cutoff}\n")
                        matchGraph.add_edge(node1, node2, homology_len = cur_homology)
   
        # while we build the initial partition give a big bonus edge for putting the homologous nodes into different partitions             
        # Adding an edge that already exists updates the edge data (in networkX graph)
        # At least one in the pair is the best similarity match for other (and second best is not too close to the best)
        # Consecutive edges are used to check but never assigned to be homologous
    for ec in matchGraph.edges():
        clear_best_match = False
        if (not G.has_edge(ec[0], ec[1])):
            for i in range (0, 2):
                best_homology = True
                best_len = matchGraph.edges[ec]['homology_len']
                for adj_node in matchGraph.neighbors(ec[i]):
                    if adj_node != ec[1 - i] and  best_len * 0.8 < matchGraph.edges[ec[i], adj_node]['homology_len']:
                        best_homology = False
                clear_best_match = clear_best_match or best_homology
            if clear_best_match:
                matchGraph.add_edge(ec[0], ec[1], weight = homology_weight, homology_len = best_len)
            else:
        #not really look like homologous node pair but still suspicious, lets just wipe the hi-c links but not prioritize splitting them to different partitions.                                
                matchGraph.add_edge(ec[0], ec[1], weight=0, homology_len = best_len)
        else:
            matchGraph.remove_edge(ec[0],ec[1])

    sys.stderr.write("Loaded match info with %d nodes and %d edges\n" % (matchGraph.number_of_nodes(), matchGraph.number_of_edges()))

    #TODO move to logging
    sys.stderr.write(f"MatchGraph edges\n")
    for d in sorted(matchGraph.edges()):
        sys.stderr.write(f"match_graph edge {d} : {matchGraph.edges[d]}  \n")

    return matchGraph

def isDiploid(matchGraph, G, component):
    hom_length = 0
    total_length = 0
    for node in component:
        total_length += G.nodes[node]['length']
        if node in matchGraph.nodes():
            for adj_node in matchGraph.neighbors(node):
                if adj_node in component and matchGraph.edges[node, adj_node]['weight'] < 0:
                    hom_length += matchGraph.edges[node, adj_node]['homology_len']
    if hom_length * 2 > total_length:
        return True
    else:
        return False