#!/usr/bin/env python

import sys
import os
import re

import fasta_util as seq

mode             = None
output_name      = None
namedict         = dict()
filenames        = []

if (len(sys.argv) > 4) and (sys.argv[1] == 'extract'):
  mode             =                 sys.argv[1]
  output_name      =                 sys.argv[2]
  namedict         = seq.readNameMap(sys.argv[3])
  filenames        =                 sys.argv[4:]

else:
  sys.stderr.write(f"usage:\n")
  sys.stderr.write(f"  {sys.argv[0]} extract <output-file> <list-of-names> <files ....>\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  writes sequences contained in <list-of-names> to stdout as uncompressed FASTA.\n")
  sys.stderr.write(f"\n")
  sys.exit(1)


outf = seq.openOutput(output_name)

for filename in filenames:
  print(f"Starting file {filename}.")
  inf  = seq.openInput(filename)

  line = inf.readline()

  while (line != ""):

    #  Read either a FASTA or FASTQ string from the input, based on the first letter
    #  of the heder line.
    #
    if   (line[0] == ">"):
      line, sName, sSeq, sQlt = seq.readFastA(inf, line)
      sName = seq.replaceName(sName, namedict, mode)

      if sName:
        outf.write(f">{sName}\n{sSeq}\n".encode())

    elif (line[0] == "@"):
      line, sName, sSeq, sQlt = seq.readFastQ(inf, line)
      sName = seq.replaceName(sName, namedict, mode)

      if sName:
        outf.write(f">{sName}\n{sSeq}\n".encode())

    else:
      print(f"Unrecognized line '{line.strip()}'", file=sys.stderr)
      line = inf.readline()

  inf.close()


if outf != None:
  outf.close()
