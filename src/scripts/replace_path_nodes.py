#!/usr/bin/env python

import sys
import re

def reverse_dir(dir):
    if dir == '>':
        return '<'
    elif dir == '<':
        return '>'

def reverse_path(path, regex_expr):
    out_path = ''
    node_names = re.findall(regex_expr, path)
    node_dirs = re.findall('[<>]', path)
    for i, node_name in (enumerate(node_names)):
        out_path = reverse_dir(node_dirs[i]) + node_name + out_path
    return out_path

def process_map_file(map_file):
    utig_map = dict()
    for line in map_file:
        line_split = line.strip().split('\t')
        utig_node = line_split[0]
        nodeID_path = line_split[1].split(':')[0]
        utig_map[utig_node] = nodeID_path
    return utig_map

in_gaf = open(sys.argv[1])
map_file = open(sys.argv[2])

utig_to_nodeIDPath = process_map_file(map_file)
for line in in_gaf:
    line_split = line.strip().split('\t')
    path = line_split[5]
    node_names = re.findall('utig1-\d+', path)
    node_dirs = re.findall('[<>]', path)
    new_path = ''
    for i in range(len(node_names)):
        node_dir = node_dirs[i]
        node_name = node_names[i]
        if node_dir == '>':
            new_path += utig_to_nodeIDPath[node_name]
        elif node_dir == '<':
            new_path += reverse_path(utig_to_nodeIDPath[node_name], '[<>]([a-zA-Z0-9-_]*)')
    line_split[5] = new_path
    print('\t'.join(line_split))
