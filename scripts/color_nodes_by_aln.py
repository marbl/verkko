#!/usr/bin/env python

import sys

# gaf from stdin
# csv to stdout

def get_color(contigs):
	h = hash(contigs)
	hue = float(hash(h) % 360)
	bright_type = h % 3
	brightness_max = 128
	brightness_min = 0
	if bright_type == 0 or bright_type == 1: brightness_max = 255
	if bright_type == 0: brightness_min = 127
	r = 0
	g = 0
	b = 0
	if hue >= 0 and hue <= 60: r = hue / 60
	if hue >= 60 and hue <= 180: r = 1
	if hue >= 180 and hue <= 240: r = (270 - hue) / 60
	if hue >= 120 and hue <= 180: g = (hue-120) / 60
	if hue >= 180 and hue <= 300: g = 1
	if hue >= 300 and hue <= 360: g = (360 - hue) / 60
	if hue >= 240 and hue <= 300: b = (240-180) / 60
	if hue >= 300 or hue <= 60: b = 1
	if hue >= 60 and hue <= 120: b = (120 - hue) / 60
	r = int(r * (brightness_max - brightness_min)) + brightness_min
	g = int(g * (brightness_max - brightness_min)) + brightness_min
	b = int(b * (brightness_max - brightness_min)) + brightness_min
	return "#" + '%02X'%r + '%02X'%g + '%02X'%b

alns_per_node = {}

for l in sys.stdin:
	parts = l.strip().split('\t')
	node = parts[0]
	contig = parts[5]
	if node not in alns_per_node: alns_per_node[node] = set()
	alns_per_node[node].add(contig)

print("node\tcontigs\tcolour")

for node in alns_per_node:
	contigs = tuple(alns_per_node[node])
	color = get_color(contigs)
	print(node + "\t" + ",".join(contigs) + "\t" + str(color))
