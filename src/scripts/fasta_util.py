#!/usr/bin/env python

import gzip
import bz2
import lzma
import re
import subprocess

#
#  Reasonably simple methods for reading FASTA and FASTQ files, munging
#  contig names, and building scaffolds.
#


##########
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
  elif ((filename[-4:] == ".bam") or
        (filename[-5:] == ".cram")):
    pipe = subprocess.Popen(["samtools", "fasta", filename], stdin=subprocess.DEVNULL, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    inf  = pipe.stdout
  else:
    inf = open(filename, mode='rt')

  return(inf)


##########
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


##########
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
  seqlines = []

  while (line) and (line[0] != '>') and (line[0] != '@'):
    seqlines.append(line.strip().upper())
    line  = inf.readline()

  sSeq = "".join(seqlines)

  return(line, sName, sSeq, sQlt)


def readFastQ(inf, line):
  sName = line[1:].strip().split()[0]
  sSeq  = inf.readline().strip().upper()
  line  = inf.readline()
  sQlt  = inf.readline().strip()
  line  = inf.readline()

  return(line, sName, sSeq, sQlt)


##########
#  Compress homopolymer runs to a single letter.
#
def homoPolyCompress(seq):
  ec  = seq[0]
  hpc = seq[0]

  for ch in seq[1:]:
    if ch != ec:
      ec   = ch
      hpc += ch

  return(hpc)


##########
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


##########
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


##########
#  Read directions for combining contigs and gaps into contigs.
#
def readScfMap(filename):
  scfmap  = dict()
  names   = dict()
  ctgname = ""
  ctglist = []

  with openInput(filename) as f:
    for l in f:
      words = re.findall(r"(\S+)", l)

      if words[0] == "path":
        ctgname = words[1]
        path    = words[2]
        ctglist = []
      elif words[0] == "end":
        scfmap[ctgname] = ctglist
        names[ctgname]  = path
      else:
        ctglist.append(words[0])

  return(scfmap, names)

