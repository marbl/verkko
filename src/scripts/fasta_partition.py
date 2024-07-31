#!/usr/bin/env python

import sys
import os
import re

import fasta_util as seq
import graph_functions as gf

mode             =  None
output_prefix    =  None
output_max_bytes = 3000000000
output_max_seqs  = 150000
min_read_length =  0
filenames        = []
output_hpc       = 0

if   (len(sys.argv) > 8) and (sys.argv[1] == 'partition'):
  mode             =             sys.argv[1]
  output_prefix    =             sys.argv[2]
  output_max_bytes =         int(sys.argv[3])
  output_max_seqs  =         int(sys.argv[4])
  min_read_length  =         int(sys.argv[5])
  output_names     = gf.str2bool(sys.argv[6])
  output_hpc       = gf.str2bool(sys.argv[7])
  filenames        =             sys.argv[8:]

else:
  sys.stderr.write(f"usage:\n")
  sys.stderr.write(f"  {sys.argv[0]} partition <output-prefix> <max-bytes> <max-seqs> <min-read-length> <files ....>\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  copies sequences from input files to output files no bigger than <max-bytes>\n")
  sys.stderr.write(f"  and containing no more than <max-seqs> sequences each.  Reads shorter than\n")
  sys.stderr.write(f"  <min-read-length> are discarded.  Output is FASTA, gzip level-1 compressed.\n")
  sys.stderr.write(f"\n")
  sys.exit(1)


outf_index       = 1
outf_bytes       = 0
outf_seqs        = 0
outf             = seq.openOutput(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz")

for filename in filenames:
  print(f"Starting file {filename}.", file=sys.stderr)
  inf  = seq.openInput(filename)

  line = inf.readline()

  while (line != ""):

    #  If the current file is big enough, close it, but delay making a new input
    #  until we have a new read to write.
    #
    if ((output_max_bytes > 0) and
        (output_max_seqs  > 0) and
        ((outf_bytes > output_max_bytes) or
         (outf_seqs  > output_max_seqs))):
      print(f"File {output_prefix}{format(outf_index, '03d')}.fasta.gz complete with {outf_bytes} bp and {outf_seqs} reads.", file=sys.stderr)
      outf_bytes  = 0
      outf_seqs   = 0
      outf_index += 1
      outf.close()
      outf        = None

    #  Read either a FASTA or FASTQ string from the input, based on the first letter
    #  of the heder line.
    #
    if   (line[0] == ">"):
      line, sName, sSeq, sQlt = seq.readFastA(inf, line)

      if output_hpc:
        sSeq = seq.homoPolyCompress(sSeq)

      if (len(sSeq) < min_read_length):
        sName = None

      if sName:
        outf_bytes += len(sSeq)
        outf_seqs  += 1

        if outf == None:
          outf        = seq.openOutput(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz")
        outf.write(f">{sName}\n{sSeq}\n".encode())
        if output_names: print(sName)

    elif (line[0] == "@"):
      line, sName, sSeq, sQlt = seq.readFastQ(inf, line)

      if output_hpc:
        sSeq = seq.homoPolyCompress(sSeq)

      if (len(sSeq) < min_read_length):
        sName = None

      if sName:
        outf_bytes += len(sSeq)
        outf_seqs  += 1

        if outf == None:
          outf        = seq.openOutput(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz")
        outf.write(f">{sName}\n{sSeq}\n".encode())
        if output_names: print(sName)

    else:
      print(f"Unrecognized line '{line.strip()}'", file=sys.stderr)
      line = inf.readline()

  inf.close()


if (output_max_bytes > 0) and (output_max_seqs  > 0):
  print(f"File {output_prefix}{format(outf_index, '03d')}.fasta.gz complete with {outf_bytes} bp and {outf_seqs} reads.", file=sys.stderr)


if outf != None:
  outf.close()
