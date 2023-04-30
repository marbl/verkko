import networkx as nx
import sys

def remove_large_tangles(G, MIN_LEN, MAX_SHORT_COMPONENT):
    shorts = set()
    for node in G.nodes():
        if G.nodes[node]['length'] < MIN_LEN:
            shorts.add(node)
    sh_G = G.subgraph(shorts)
    nodes_deleted = 0
    components_deleted = 0
    to_delete = []
    for comp in nx.connected_components(sh_G):
        if len(comp) > MAX_SHORT_COMPONENT:
            components_deleted += 1
            for e in comp:
                nodes_deleted += 1
                to_delete.append(e)
    G.remove_nodes_from(to_delete)
    sys.stderr.write(f'Removed {components_deleted} short nodes components and {nodes_deleted} short nodes. New '
                     f'number of nodes {G.number_of_nodes()}\n')
    return set(to_delete)

def load_indirect_graph(gfa_file, G):
    translate = open(gfa_file, 'r')
    for line in translate:
        if "#" in line:
            continue
        line = line.strip().split()

        if line[0] == "S":
            G.add_node(line[1], length=int(line[3][5:]), coverage=float(line[5][5:]))
        elif line[0] == "L":
            if line[1] not in G or line[3] not in G:
                sys.stderr.write("Warning, skip link between nodes not in graph:%s" % (line))
                sys.exit(1)
            G.add_edge(line[1], line[3])
