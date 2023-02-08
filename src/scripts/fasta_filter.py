#!/usr/bin/env python

import sys
import os
import re

import fasta_util as seq

mode             =  None
output_prefix    =  None
output_max_bytes = 3000000000
output_max_seqs  = 150000
min_read_length =  0
filenames        = []
output_hpc       = 0

if   (len(sys.argv) > 6) and (sys.argv[1] == 'filter'):
  mode             =     sys.argv[1]
  output_prefix    =     sys.argv[2]
  output_max_bytes = int(sys.argv[3])
  output_max_seqs  = int(sys.argv[4])
  min_read_length =  int(sys.argv[5])
  filenames        =     sys.argv[6:]
  output_hpc       = 0

else:
  sys.stderr.write(f"usage:\n")
  sys.stderr.write(f"  {sys.argv[0]} partition <output-prefix> <max-bytes> <max-seqs> <min-read-length> <files ....>\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  copies sequences from input files to output files no bigger than <max-bytes>\n")
  sys.stderr.write(f"  and containing no more than <max-seqs> sequences each.  Reads shorter than\n")
  sys.stderr.write(f"  <min-read-length> are discarded.  Output is FASTA, gzip level-1 compressed.\n")
  sys.stderr.write(f"\n")
  sys.exit(1)


outf = seq.openOutput(output_name)
