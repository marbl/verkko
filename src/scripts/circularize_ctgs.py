#!/usr/bin/env python
import sys
import argparse
import parasail
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser(description="Circularize contigs")
parser.add_argument("contigs", nargs='?', help="Input FASTA file (by default reading FASTA from stdin)")
parser.add_argument("-p", "--penalty", type=int, default=100, help="Mismatch/indel penalty (default: 100bp)")
parser.add_argument("-f", "--flank", type=float, default=0.60, help="Maximum fraction allowed to be trimmed (default: 0.60)")
parser.add_argument("-o", "--output", help="Output FASTA (by default prints to stdout)")
parser.add_argument("--min-ovl", type=int, default=0, help="Minimal length of overlap (default 0)")
args = parser.parse_args()

if args.contigs:
    print("Reading from file", args.contigs, file=sys.stderr)
    i_handle = open(args.contigs, 'r')
else:
    print("Reading from stdin", file=sys.stderr)
    i_handle = sys.stdin

penalty = args.penalty
flank = args.flank

print('Mismatch/indel penalty', penalty, file=sys.stderr)
print('Using Parasail', file=sys.stderr)

def parasail_ends(s1, s2):
    m = parasail.matrix_create("ACGT", 1, -penalty)
    res = parasail.sg_qb_de(s1, s2, penalty, penalty, m)
    return res.end_ref + 1, res.score

def overlap_and_trim(s1):
    f1 = min(len(s1), flank)
    f2 = min(len(s1), flank)
    print('Aligning', file=sys.stderr)
    pos, score = parasail_ends(str(s1[-f1:]), str(s1[:f2]))

    if score > 0:
        if pos < args.min_ovl:
            print('Overlap', pos, 'was smaller than threshold', args.min_ovl, file=sys.stderr)
            return s1
        else:
            print('Overlap', pos, file=sys.stderr)
            print('Score', score, file=sys.stderr)

            return s1[pos:]  
    else:
        print('Overlap score', score, 'was below zero', file=sys.stderr) 
        return s1

def make_sequence(seq):
    return overlap_and_trim(seq)

contig_dict = dict()
for r in SeqIO.parse(i_handle, "fasta"):
    contig_dict[r.name] = r.seq

if args.contigs:
    i_handle.close()

records = []
for n in contig_dict:
    print('=======================================', file=sys.stderr)
    flank = int(len(contig_dict[n]) * flank)
    print('Processing', n, 'with flank set to', flank, 'and length', len(contig_dict[n]), file=sys.stderr)
    seq = make_sequence(contig_dict[n])
    records.append(SeqRecord(seq, id=n, description=''))

if args.output:
    print("Writing to file", args.output, file=sys.stderr)
    o_handle = open(args.output, 'w')
else:
    print("Writing to stdin", file=sys.stderr)
    o_handle = sys.stdout

SeqIO.write(records, o_handle, 'fasta')

if args.output:
    o_handle.close()
