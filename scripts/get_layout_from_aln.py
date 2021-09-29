#!/usr/bin/python

import sys
import random

ref_file = sys.argv[1]
paf_file = sys.argv[2]
read_names_file = sys.argv[3]
read_files = sys.argv[4:]

def add_lines(read_name, current_lines, lines_per_contig, read_order):
	if len(current_lines) == 0: return
	if read_name not in read_order: return
	max_mapq_lines = []
	max_mapq = 0
	for l in current_lines:
		parts = l.split('\t')
		if float(int(parts[3]) - int(parts[2])) / float(int(parts[1])) < 0.98: continue
		if float(int(parts[8]) - int(parts[7])) / float(int(parts[1])) > 1.02: continue

		if int(parts[11]) < max_mapq: continue
		if int(parts[11]) > max_mapq:
			max_mapq_lines = []
			max_mapq = int(parts[11])
		if int(parts[11]) == max_mapq: max_mapq_lines.append(l)
	if len(max_mapq_lines) == 0: return
	assert len(max_mapq_lines) >= 1
	# if ambiguously aligned, randomly pick a position
	# random_index = random.randint(0, len(max_mapq_lines))
	random_index = 0
	l = max_mapq_lines[random_index]
	parts = l.split('\t')
	left_clip = int(parts[2])
	right_clip = int(parts[1]) - int(parts[3])
	start_pos = int(parts[7])
	end_pos = int(parts[8])
	if parts[4] == "-":
		end_pos += left_clip
		start_pos -= right_clip
		(start_pos, end_pos) = (end_pos, start_pos)
	else:
		start_pos -= left_clip
		end_pos += right_clip
	lines_per_contig[parts[5]].append((min(int(start_pos), int(end_pos)), "read\t" + str(read_order[read_name]) + "\tanchor\t0\thang\t0\t0\tposition", start_pos, end_pos))

read_order = {}
next_id = 1
for read_file in read_files:
	read_name = ""
	read_len = 0
	with open(read_file) as f:
		for l in f:
			if l[0] == '>':
				if read_len >= 1000 and read_name != "":
					read_order[read_name] = next_id
					next_id += 1
				read_name = l[1:].split(' ')[0].strip()
				read_len = 0
			else:
				read_len += len(l.strip())
	if read_len >= 1000 and read_name != "":
		read_order[read_name] = next_id
		next_id += 1

with open(read_names_file, "w") as f:
	for name in read_order: f.write(str(name) + "\t" + str(read_order[name]) + "\n")

lines_per_contig = {}
contig_len = {}
contig_name = ""
current_len = 0
with open(ref_file) as f:
	for l in f:
		if l[0] == '>':
			if contig_name != "":
				lines_per_contig[contig_name] = []
				contig_len[contig_name] = current_len
			contig_name = l[1:].strip()
			current_len = 0
		else:
			current_len += len(l.strip())
if contig_name != "":
	lines_per_contig[contig_name] = []
	contig_len[contig_name] = current_len

current_name = ""
current_lines = []
with open(paf_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if parts[0] != current_name:
			add_lines(current_name, current_lines, lines_per_contig, read_order)
			current_lines = []
			current_name = parts[0]
		current_lines.append(l.strip())

add_lines(current_name, current_lines, lines_per_contig, read_order)
current_lines = []
current_name = parts[0]
contig_ids = {}
next_id = 0

for contig in lines_per_contig:
	if len(lines_per_contig[contig]) == 0: continue
	lines_per_contig[contig].sort(key=lambda x: x[0])
	start_pos = contig_len[contig]
	end_pos = 0
	for line in lines_per_contig[contig]:
		start_pos = min(start_pos, line[2])
		start_pos = min(start_pos, line[3])
		end_pos = max(end_pos, line[2])
		end_pos = max(end_pos, line[3])
	contig_ids[contig] = next_id
	next_id += 1
	print("tig\t-1")
	print("len\t" + str(end_pos - start_pos))
	print("cns")
	print("qlt")
	print("trimBgn\t0")
	print("trimEnd\t" + str(end_pos - start_pos))
	print("suggestRepeat\tF")
	print("suggestBubble\tF")
	print("suggestCircular\tF")
	print("circularLength\t0")
	print("numChildren\t" + str(len(lines_per_contig[contig])))
	for line in lines_per_contig[contig]:
		print(line[1] + "\t" + str(line[2] - start_pos) + "\t" + str(line[3] - start_pos))
	print("tigend")

with open("tig_names.txt", "w") as f:
   for name in contig_ids: f.write(str(contig_ids[name]) + "\t" + str(name) + "\n") 
