#!/usr/bin/env python

import sys
import os
import re

import fasta_util as seq

mode             = None
output_name      = None
namedict         = dict()
scfmap           = []
pieces           = dict()
filenames        = []

if   (len(sys.argv) > 4) and (sys.argv[1] == 'rename'):
  mode             =                 sys.argv[1]
  output_name      =                 sys.argv[2]
  namedict         = seq.readNameMap(sys.argv[3])
  filenames        =                 sys.argv[4:]

elif (len(sys.argv) > 5) and (sys.argv[1] == 'combine'):
  mode             =                 sys.argv[1]
  output_name      =                 sys.argv[2]
  namedict         = seq.readNameMap(sys.argv[3])
  scfmap,_         =  seq.readScfMap(sys.argv[4])
  filenames        =                 sys.argv[5:]

else:
  sys.stderr.write(f"usage:\n")
  sys.stderr.write(f"  {sys.argv[0]} rename    <output-file>   <name-map>                <files ...>\n")
  sys.stderr.write(f"  {sys.argv[0]} combine   <output-file>   <name-map> <scaffold-map> <files ...>\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'rename' changes the name of sequences according to 'name-map' or leaves\n")
  sys.stderr.write(f"  alone if not in the map and writes the output to stdout.  Output is uncompressed FASTA.\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'combine' performs 'rename' and then combines sequences to scaffolds as\n")
  sys.stderr.write(f"  described in the <scaffold-map>.  Output is uncompressed FASTA.\n")
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

      if scfmap:
        pieces[sName] = sSeq
      elif sName:
        outf.write(f">{sName}\n{sSeq}\n".encode())

    elif (line[0] == "@"):
      line, sName, sSeq, sQlt = seq.readFastQ(inf, line)
      sName = seq.replaceName(sName, namedict, mode)

      if scfmap:
        pieces[sName] = sSeq
      elif sName:
        outf.write(f">{sName}\n{sSeq}\n".encode())

    else:
      print(f"Unrecognized line '{line.strip()}'", file=sys.stderr)
      line = inf.readline()

  inf.close()


if scfmap:  #  For 'combine' combine the contig pieces and gaps into scaffolds.
  print(f"\nWriting output.")

  for clist in scfmap:
    seq = ""
    prev = ""
    for piece in scfmap[clist]:
      numn = re.match(r"\[N(\d+)N]", piece)
      if numn:
        if not seq:
           print(f"ERROR:piece {prev} missing from gapped contig {clist}.", file=sys.stderr)
           sys.exit(1)
        seq += "N" * int(numn[1])
      elif piece in pieces:
        seq += pieces[piece]
      elif seq:
        print(f"ERROR: piece {piece} missing from gapped contig {clist}.", file=sys.stderr)
        sys.exit(1)
      prev = piece

    if seq:
      outf.write(f">{clist}\n".encode())
      outf.write(f"{seq}\n".encode())
    else:
      print(f"ERROR: contig {clist} is empty.", file=sys.stderr)

  print(f"")
  print(f"Finished.")


if outf != None:
  outf.close()



