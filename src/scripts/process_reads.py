#!/usr/bin/env python

import sys
import os
import gzip
import bz2
import lzma
import re



#  Given the next line in the file, read a FASTA/FASTQ formatted sequence and
#  return the ident line, the sequence and quality values (or empty string if
#  FASTA).
#
#  Also returns the line after the sequence.
#
def readFastA(inf, line):
  sName = line[1:].strip().split()[0]
  sSeq = ""
  sQlt = ""

  line  = inf.readline()

  while (line) and (line[0] != '>') and (line[0] != '@'):
    sSeq += line.strip().upper()
    line  = inf.readline()

  return(line, sName, sSeq, sQlt)


def readFastQ(inf, line):
  sName = line[1:].strip().split()[0]
  sSeq  = inf.readline().strip().upper()
  line  = inf.readline()
  sQlt  = inf.readline().strip()
  line  = inf.readline()

  return(line, sName, sSeq, sQlt)


def homoPolyCompress(seq):
  ec  = seq[0]
  hpc = seq[0]

  for ch in seq[1:]:
    if ch != ec:
      ec   = ch
      hpc += ch

  return(hpc)


#  Open an input.  gzip/xz/bzip2/stdin/uncompressed all supported.  (But not
#  compressed stdin, since we use the filename to decide how to open it.)
#
def openInput(filename):
  if   (filename      == "-"):
    inf = sys.stdin
  elif (filename[-3:] == ".gz"):
    inf = gzip.open(filename, mode='rt')
  elif (filename[-3:] == ".xz"):
    inf = lzma.open(filename, mode='rt')
  elif (filename[-4:] == ".bz2"):
    inf = bz2.open(filename, mode='rt')
  else:
    inf = open(filename, mode='rt')

  return(inf)


#  Open an output.  gzip/xz/bzip2/stdin/uncompressed all supported.  (But not
#  compressed stdout, since we use the filename to decide how to open it.)
#
#  stdout is reopened as binary, for compatibility with the compressed
#  writers.
#
def openOutput(filename):
  if   (filename      == "-"):
    #otf = sys.stdout
    otf = os.fdopen(sys.stdout.fileno(), 'wb', closefd=False)
  elif (filename[-3:] == ".gz"):
    otf = gzip.open(filename, mode='wb', compresslevel=1)
  elif (filename[-3:] == ".xz"):
    otf = lzma.open(filename, mode='wb', compresslevel=1)
  elif (filename[-4:] == ".bz2"):
    otf = bz2.open(filename, mode='wb', compresslevel=1)
  else:
    otf = open(filename, mode='wb')

  return(otf)


#  Read a tig ID to tig name map.
#
def readNameMap(filename):
  namedict = dict()

  with openInput(filename) as f:
    for l in f:
      w = l.strip().split()
      if len(w) == 1:
        namedict[w[0]] = w[0]
      else:
        namedict[w[0]] = w[1]

  return(namedict)


#  Replace the 'name' (excluding any 'tig0..' prefix) with whatever is in the
#  dictionary.  If nothing matches, the original name is retained.
#
def replaceName(name, namedict, mode):
  if   mode == 'partition':
    return(name)
  elif mode == 'extract':
    return(namedict.get(name))
  elif mode == 'rename' or mode == 'combine':
    return(namedict.get(re.sub('tig0+', '', name), name))
  else:
    return(name)

#  Return true if 'name' is in the dictionary.
#
def existsInDict(name, namedict):
  return(name in namedict)


def readScfMap(filename):
  scfmap  = dict()
  ctgname = ""
  ctglist = []

  with openInput(filename) as f:
    for l in f:
      words = re.findall(r"(\S+)", l)

      if words[0] == "path":
        ctgname = words[1]
        ctglist = []
        print(f"{ctgname}")
      elif words[0] == "end":
        scfmap[ctgname] = ctglist
        print(f"")
      else:
        ctglist.append(words[0])
        print(f" - {words[0]}")

  return(scfmap)


mode             = None
filenames        = []

namedict         = dict()
scfmap           = []
pieces           = dict()

output_prefix    = None
output_max_bytes = 0
output_max_seqs  = 0

output_hpc       = 0

outf_index       = 1
outf_bytes       = 0
outf_seqs        = 0
outf             = None

min_read_length  = 0

if   (len(sys.argv) > 6) and (sys.argv[1] == 'partition'):
  mode             =     sys.argv[1]
  output_prefix    =     sys.argv[2]
  output_max_bytes = int(sys.argv[3])
  output_max_seqs  = int(sys.argv[4])
  min_read_length =  int(sys.argv[5])
  filenames        =     sys.argv[6:]
  output_hpc       = 0

elif (len(sys.argv) > 4) and (sys.argv[1] == 'extract'):
  mode             =             sys.argv[1]
  output_name      =             sys.argv[2]
  namedict         = readNameMap(sys.argv[3])
  filenames        =             sys.argv[4:]
  output_hpc       = 0

elif (len(sys.argv) > 4) and (sys.argv[1] == 'rename'):
  mode             =             sys.argv[1]
  output_name      =             sys.argv[2]
  namedict         = readNameMap(sys.argv[3])
  filenames        =             sys.argv[4:]

elif (len(sys.argv) > 5) and (sys.argv[1] == 'combine'):
  mode             =             sys.argv[1]
  output_name      =             sys.argv[2]
  namedict         = readNameMap(sys.argv[3])
  scfmap           =  readScfMap(sys.argv[4])
  filenames        =             sys.argv[5:]

else:
  sys.stderr.write(f"usage:\n")
  sys.stderr.write(f"  {sys.argv[0]} partition <output-prefix> <max-bytes> <max-seqs>    <files ....>\n")
  sys.stderr.write(f"  {sys.argv[0]} extract   <output-file>   <list-of-names>           <files ....>\n")
  sys.stderr.write(f"  {sys.argv[0]} rename    <output-file>   <name-map>                <files ...>\n")
  sys.stderr.write(f"  {sys.argv[0]} combine   <output-file>   <name-map> <scaffold-map> <files ...>\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'partition' creates output files no bigger than <max-bytes> and containing\n")
  sys.stderr.write(f"   no more than <max-seqs> sequences each.  Output is FASTA, gzip level-1\n")
  sys.stderr.write(f"   compressed.\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'extract' writes sequences contained in <list-of-names> to stdout as uncompressed.\n")
  sys.stderr.write(f"  FASTA\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'rename' changes the name of sequences according to 'name-map' or leaves\n")
  sys.stderr.write(f"  alone if not in the map and writes the output to stdout.  Output is uncompressed FASTA.\n")
  sys.stderr.write(f"\n")
  sys.stderr.write(f"  'combine' performs 'rename' and then combines sequences to scaffolds as\n")
  sys.stderr.write(f"  described in the <scaffold-map>.  Output is uncompressed FASTA.\n")
  sys.stderr.write(f"\n")
  sys.exit(1)


if output_prefix:
  outf = openOutput(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz")
else:
  outf = openOutput(output_name)


for filename in filenames:
  print(f"Starting file {filename}.")
  inf  = openInput(filename)

  line = inf.readline()

  while (line != ""):

    #  If the current file is big enough, close it, but delay making a new input
    #  until we have a new read to write.
    #
    if ((output_prefix) and
        (output_max_bytes > 0) and
        (output_max_seqs  > 0) and
        ((outf_bytes > output_max_bytes) or
         (outf_seqs  > output_max_seqs))):
      print(f"File {output_prefix}{format(outf_index, '03d')}.fasta.gz complete with {outf_bytes} bp and {outf_seqs} reads.")
      outf_bytes  = 0
      outf_seqs   = 0
      outf_index += 1
      outf.close()
      outf        = None

    #  Read either a FASTA or FASTQ string from the input, based on the first letter
    #  of the heder line.
    #
    if   (line[0] == ">"):
      line, sName, sSeq, sQlt = readFastA(inf, line)
      sName = replaceName(sName, namedict, mode)

      if output_hpc:
        sSeq = homoPolyCompress(sSeq)

      if (mode =='partition') and (len(sSeq) < min_read_length):
        sName = None

      if scfmap:
        print(f"SAVE {sName}")
        pieces[sName] = sSeq
      elif sName:
        #print(f"Write '{sName}' of length {len(sSeq)} to index {format(outf_index, '03d')}")
        outf_bytes += len(sSeq)
        outf_seqs  += 1

        if outf == None:
          outf        = gzip.open(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz", mode="wb", compresslevel=1)
        outf.write(f">{sName}\n{sSeq}\n".encode())

    elif (line[0] == "@"):
      line, sName, sSeq, sQlt = readFastQ(inf, line)
      sName = replaceName(sName, namedict, mode)

      if output_hpc:
        sSeq = homoPolyCompress(sSeq)

      if (mode =='partition') and (len(sSeq) < min_read_length):
        sName = None

      if scfmap:
        print(f"SAVE {sName}")
        pieces[sName] = sSeq
      elif sName:
        #print(f"Write '{sName}' of length {len(sSeq)} to index {format(outf_index, '03d')}")
        outf_bytes += len(sSeq)
        outf_seqs  += 1

        if outf == None:
          outf        = gzip.open(f"{output_prefix}{format(outf_index, '03d')}.fasta.gz", mode="wb", compresslevel=1)
        outf.write(f">{sName}\n{sSeq}\n".encode())

    else:
      print(f"Unrecognized line '{line.strip()}'")
      line = inf.readline()

  inf.close()


if scfmap:
  print(f"\nWriting output.")
  print(f"")

  for clist in scfmap:
    print(f">{clist}")

    seq = ""
    for piece in scfmap[clist]:
      numn = re.match(r"\[N(\d+)N]", piece)

      if numn:
        print(f"  - {numn[1]} Ns")
        seq += "N" * int(numn[1])
      else:
        print(f"  - {piece}")
        seq += pieces[piece]

    outf.write(f">{clist}\n".encode())
    outf.write(f"{seq}\n".encode())

  print(f"")
  print(f"Finished.")


if ((output_prefix) and
    (output_max_bytes > 0) and
    (output_max_seqs  > 0)):
  print(f"File {output_prefix}{format(outf_index, '03d')}.fasta.gz complete with {outf_bytes} bp and {outf_seqs} reads.")


if outf != None:
  outf.close()
