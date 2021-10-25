#!/usr/bin/env python

import sys

pat_kmer_file = sys.argv[1]
mat_kmer_file = sys.argv[2]
# graph from stdin
# csv to stdout

def revcomp(s):
	comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
	return "".join(comp[c] for c in s[::-1])

def kmer_to_int(s):
	result = 0
	for i in range(0, len(s)):
		result <<= 2
		if s[i] == 'A':
			result += 0
		elif s[i] == 'C':
			result += 1
		elif s[i] == 'G':
			result += 2
		elif s[i] == 'T':
			result += 3
		elif s[i] == "N":
			return 0
		else:
			assert False
	return result

pat_kmers = set()
mat_kmers = set()

kmer_size = None

sys.stderr.write("read paternal k-mers\n")

with open(pat_kmer_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if kmer_size is None: kmer_size = len(parts[0])
		assert len(parts[0]) == kmer_size
		pat_kmers.add(kmer_to_int(parts[0]))

sys.stderr.write("read maternal k-mers\n")

with open(mat_kmer_file) as f:
	for l in f:
		parts = l.strip().split('\t')
		if kmer_size is None: kmer_size = len(parts[0])
		assert len(parts[0]) == kmer_size
		mat_kmers.add(kmer_to_int(parts[0]))

assert kmer_size is not None

sys.stderr.write("read nodes\n")

print("node,pat,mat,merge")

for l in sys.stdin:
	parts = l.strip().split('\t')
	if parts[0] != 'S': continue
	mat_count = 0
	pat_count = 0
	fw_kmer = 0
	seqs = parts[2].split('N')
	for seq in seqs:
		if len(seq) < kmer_size: continue
		for i in range(0, kmer_size):
			fw_kmer <<= 2
			if seq[i] == 'A':
				fw_kmer += 0
			elif seq[i] == 'C':
				fw_kmer += 1
			elif seq[i] == 'G':
				fw_kmer += 2
			elif seq[i] == 'T':
				fw_kmer += 3
			else:
				assert False
		if fw_kmer in mat_kmers: mat_count += 1
		if fw_kmer in pat_kmers: pat_count += 1
		for i in range(kmer_size, len(seq)):
			fw_kmer <<= 2
			if seq[i] == 'A':
				fw_kmer += 0
			elif seq[i] == 'C':
				fw_kmer += 1
			elif seq[i] == 'G':
				fw_kmer += 2
			elif seq[i] == 'T':
				fw_kmer += 3
			else:
				assert False
			fw_kmer = fw_kmer % (pow(4, kmer_size))
			if fw_kmer in mat_kmers: mat_count += 1
			if fw_kmer in pat_kmers: pat_count += 1
	seqs = revcomp(parts[2]).split('N')
	for seq in seqs:
		if len(seq) < kmer_size: continue
		for i in range(0, kmer_size):
			fw_kmer <<= 2
			if seq[i] == 'A':
				fw_kmer += 0
			elif seq[i] == 'C':
				fw_kmer += 1
			elif seq[i] == 'G':
				fw_kmer += 2
			elif seq[i] == 'T':
				fw_kmer += 3
			else:
				assert False
		if fw_kmer in mat_kmers: mat_count += 1
		if fw_kmer in pat_kmers: pat_count += 1
		for i in range(kmer_size, len(seq)):
			fw_kmer <<= 2
			if seq[i] == 'A':
				fw_kmer += 0
			elif seq[i] == 'C':
				fw_kmer += 1
			elif seq[i] == 'G':
				fw_kmer += 2
			elif seq[i] == 'T':
				fw_kmer += 3
			else:
				assert False
			fw_kmer = fw_kmer % (pow(4, kmer_size))
			if fw_kmer in mat_kmers: mat_count += 1
			if fw_kmer in pat_kmers: pat_count += 1
	print(parts[1] + "," + str(pat_count) + "," + str(mat_count) + "," + str(pat_count) + "-" + str(mat_count))

