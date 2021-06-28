#!/usr/bin/python

import sys

bridge_node_file = sys.argv[1]

def canon(path):
    pathstr = "".join(path)
    reverse = [('>' if n[0] == '<' else '<') + n[1:] for n in path[::-1]]
    revstr = "".join(reverse)
    if pathstr < revstr: return path
    return reverse

bridge_nodes = set()

with open(bridge_node_file) as f:
    for l in f:
        nodes = l.strip().split(',')
        for node in nodes:
            if node.strip() == "": continue
            bridge_nodes.add(node.strip())

paths = {}

for l in sys.stdin:
    pathstr = l.strip().split('\t')[0] + '>'
    last_sep = 0
    count = 0
    bridge = []
    for i in range(0, len(pathstr)):
        if pathstr[i] == '>' or pathstr[i] == '<':
            node = pathstr[last_sep+1:i]
            if node in bridge_nodes:
                count += 1
            if count >= 1: bridge.append(pathstr[last_sep:i])
            if node in bridge_nodes and count >= 2:
                bridge = canon(bridge)
                print("".join(bridge))
                connection = tuple(canon([bridge[0], bridge[-1]]))
                if connection not in paths: paths[connection] = 0
                paths[connection] += 1
                bridge = [pathstr[last_sep:i]]
            last_sep = i

print("")
paths = [(paths[p], p) for p in paths]
paths.sort(key=lambda x: -x[0])
for p in paths:
    print(str(p[0]) + "\t" + str(p[1]))
