#!/usr/bin/env python
import sys
import re
from Bio import SeqIO
import gzip
from mimetypes import guess_type
from functools import partial

def name(r):
    return re.sub('tig0+', '', r.id)

def change_name(r):
    #print("Changing name from %s to %s" %(name(r), mapping[r.name]))
    return SeqIO.SeqRecord(r.seq, id=mapping[name(r)], description="")

def gz_open(fn, mode):
    encoding = guess_type(fn)[1]  # uses file extension
    if encoding is None:
        print("Working with text file")
        return open(fn, mode=mode)
    elif encoding == 'gzip':
        print("Working with gzipped file")
        if mode == 'r':
            mode = 'rt'
        elif mode == 'w':
            mode = 'wt'
        else:
            raise ValueError('Unknown mode "{}"'.format(mode))
        return gzip.open(fn, mode=mode)
    else:
        raise ValueError('Unknown file encoding: "{}"'.format(encoding))

def read_map(fn):
    d = dict()
    with open(fn, "r") as f:
        for l in f:
            split = l.split()
            d[split[0]] = split[1]
    return d

if len(sys.argv) < 4:
    print("Usage: %s <fasta file> <name map ('from to')> <output>" % sys.argv[0])
    sys.exit(1)

mapping = read_map(sys.argv[2])
if sys.argv[1] == "-":
   fin=sys.stdin
else:
   fin=sys.argv[1]
 
with gz_open(fin, 'r') as i_handle:
    input_seq_iterator = SeqIO.parse(i_handle, "fasta")
    with gz_open(sys.argv[3], 'w') as o_handle:
        SeqIO.write((change_name(record) for record in input_seq_iterator if name(record) in mapping), o_handle, "fasta")
