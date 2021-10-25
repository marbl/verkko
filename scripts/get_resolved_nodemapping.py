#!/usr/bin/env python

import fileinput

# gfa from stdin
# mapping to stdout

for l in fileinput.input():
	if l[0] == "S":
		name = l.strip().split('\t')[1]
		nameparts = name.split('_')[2:]
		namestr = []
		assert len(nameparts) >= 1
		for part in nameparts:
			partname = (">" if part[-1] == "f" else "<") + part[:-2]
			namestr.append(partname)
		print(name + "\t" + "".join(namestr))
