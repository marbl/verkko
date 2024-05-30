
import sys
import re
import shutil
import os
import networkx as nx
from numpy import argmax
import graph_functions
import logging


#All functions expect directed nodes in utig4-234+ format, paths may contain gapped node

#path nodes stored with orientation, "utig4-234+" tsv format
class PathStorage:
    def __init__(self):        
        #from IDs to array of oriented+- nodes
        self.paths = {}
        #ignoring overlaps but who cares
        self.path_lengths = {}
        self.hap_labels = {}

    def getPathById(self, path_id):
        return self.paths[path_id]
    
    def getLength(self, path_id):
        return self.path_lengths[path_id]

    def getPathIds(self):
        return self.paths.keys()

#should not be used
    def getEdgeSequenceById(self, path_id):
        res = []
        for edge in self.paths[path_id]:
            #removing +- and gaps
            if edge[0] != "N":
                res.append(edge[:-1])
        return res
    
    def getPathString(self, path_id):
        return ",".join(self.paths[path_id])
    
    def getPathTsv(self, path_id):
        return path_id + "\t" + ",".join(self.paths[path_id]) + "\t" + self.hap_labels[path_id]
    
    def getPathGaf(self, path_id):
        return path_id + "\t" + graph_functions.tsv2gaf(",".join(self.paths[path_id])) + "\t" + self.hap_labels[path_id]

    def addPath(self, line, G):
#        print (line)
        arr = line.strip().split()
        if len(arr) < 3:
            print (f"wrong path {line}")
            exit()
        separators = ">|<|,"
        edges = re.split(separators, arr[1])
        self.hap_labels[arr[0]] = arr[2]
        total_l = 0
        for edge in edges:
            node = edge
            if node in G.nodes:
                total_l += G.nodes[node]['length']
        edges = list(filter(None, edges))
#        print (edges)
        self.paths[arr[0]] = edges
        self.path_lengths[arr[0]] = total_l

    def addPathWithId(self, id, path, G):
        total_l = 0
        for edge in path:
            node = edge
            if node in G.nodes:
                total_l += G.nodes[node]['length']
        self.paths[id] = path
        self.path_lengths[id] = total_l
    
    def getLabel(self, path_id):
        return self.hap_labels[path_id]
