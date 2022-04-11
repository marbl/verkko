#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
from statistics import median_low

def segment_length(split_line):
    if split_line[2] != '*':
        return len(split_line[2])
    else:
        for s in split_line:
            if s.startswith("LN:i:"):
                return int(s.split(":")[2])
        assert(False)

parser = argparse.ArgumentParser(description="Integrate node sequences into GFA")
parser.add_argument("coverage", help="File with average segment coverage values")
parser.add_argument("gfa", nargs='?', help="GFA file specifying graph structure (by default reading from stdin)")
parser.add_argument("--allow-absent", action="store_true")
args = parser.parse_args()

seg_cov=dict()
with open(args.coverage, 'r') as f:
    for l in f:
        s = l.split()
        if s[1] == "coverage":
            continue
        seg_cov[s[0]] = float(s[1])

if args.gfa:
    print("Reading graph from", args.gfa, file=sys.stderr)
    stream = open(args.gfa, 'r')
else:
    print("Reading graph from stdin", file=sys.stderr)
    stream = sys.stdin

for l in stream:
    if l.startswith('S\t'):
        s = l.split()
        seg = s[1]
        s = [item for item in s if not re.match("^(RC:i:|FC:i:|ll:f:).*", item)]
        length = segment_length(s)
        assert seg in seg_cov or args.allow_absent
        cov = seg_cov[seg] if seg in seg_cov else 0.
        print('%s\tRC:i:%d\tll:f:%.1f' % ('\t'.join(s), int(length * cov), cov))
    else:
        print(l.strip())

if args.gfa:
    stream.close()
