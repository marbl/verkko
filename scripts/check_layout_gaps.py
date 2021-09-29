#!/usr/bin/python

import fileinput

# layout from stdin
# gaps to stdout

def check_alns(current_tig, current_len, current_alns):
	assert current_tig is not None
	assert current_len > 0
	current_alns.sort(key=lambda x: x[0])
	last_end = 0
	for i in range(0, len(current_alns)):
		current_start = current_alns[i][0]
		current_end = current_alns[i][1]
		if current_start > last_end:
			print(current_tig + "\t" + str(last_end) + "\t" + str(current_start) + " (len " + str(current_len) + ")")
		if current_start < 0 or current_end > current_len:
			print("layout outside contig: " + current_tig + "\t" + current_start + "\t" + current_end + "\t" + " (len " + str(current_len) + ")")
		last_end = max(last_end, current_end)
	if last_end < current_len:
		print(current_tig + "\t" + str(last_end) + "\t" + str(current_len) + " (len " + str(current_len) + ")")

current_tig = None
current_len = 0
current_alns = []

for l in fileinput.input():
	parts = l.strip().split('\t')
	if len(parts) == 0: continue
	if parts[0] == "tig":
		if current_len > 0: check_alns(current_tig, current_len, current_alns)
		current_tig = parts[1]
		current_len = 0
		current_alns = []
	if parts[0] == "len":
		assert current_len == 0
		assert current_tig is not None
		assert len(current_alns) == 0
		current_len = int(parts[1])
	if parts[0] == "read":
		start = int(parts[8])
		end = int(parts[9])
		if end < start: (start, end) = (end, start)
		current_alns.append((start, end))

if current_len > 0: check_alns(current_tig, current_len, current_alns)
