#!/usr/bin/python

import sys
import re
orient_node = sys.argv[1]

for line in sys.stdin:
    l = line.strip() + '>'
    if len(l) == 0: continue
    if orient_node not in l: continue
    assert l[0] == '>' or l[0] == '<'
    parts = []
    indices = []
    last_sep = 0
    for i in range(1, len(l)):
        if l[i] != '>' and l[i] != '<': continue
        parts.append(l[last_sep:i])
        last_sep = i
        if parts[-1][1:] == orient_node:
            indices.append(len(parts)-1)
    for index in indices:
        if index == 0 or index == len(parts)-1: continue
        triplet = parts[index-1:index+2]
        if triplet[1][0] == '<':
            triplet = triplet[::-1]
            for i in range(0, len(triplet)):
                triplet[i] = ('>' + triplet[i][1:]) if (triplet[i][0] == '<') else ('<' + triplet[i][1:])
        print(''.join(triplet))

