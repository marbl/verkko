#!/usr/bin/env python

import sys

readnamefile = sys.argv[1]

readnames = set()

with open(readnamefile) as f:
	for l in f:
		readnames.add(l.strip().split('\t')[0].split(' ')[0])

printing = False

for l in sys.stdin:
    if l[0] == '>':
        name = l.strip()[1:].split(" ")[0]
        valid_name = name in readnames
        printing = valid_name
    if printing:
	print(l.strip())
