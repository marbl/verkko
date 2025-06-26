#!/usr/bin/env python
import pysam
import sys
import fasta_util as seq
import re

# rename a bam based on the input scfmap and tig name mapping files
# in the case of scaffolds, offset the start coordinate by the start of the sequence in the scaffold
#
prefixes = [ 'hifi_', 'ont_' ]
bamfile  = pysam.AlignmentFile("-", "rb")
layout   = sys.argv[1]
gap_info = sys.argv[2]
scfmap   = seq.readScfMap(sys.argv[3])
namedict = seq.readNameMap(sys.argv[4])
lens     = dict()
offsets  = dict()  
names    = dict()
tigtohap = dict()
readtorg = dict()

with open(layout) as f:
   for l in f:
      parts = l.strip().split('\t')
      if len(parts) > 2:
         readtorg[parts[0]] = "LA" if (int(parts[-1]) == 0 or len(parts) == 3) else "UL"
with open(gap_info) as f:
   for l in f:
      parts = l.strip().split('\t')
      if parts[0] in readtorg and readtorg[parts[0]] == "LA":
          readtorg[parts[0]] = "UL-gap"

for filename in sys.argv[5:]:
  sys.stderr.write("Starting file %s\n"%(filename))
  inf  = seq.openInput(filename)

  line = inf.readline()

  while (line != ""):
     if   (line[0] == ">"):
        line, sName, sSeq, sQlt = seq.readFastA(inf, line)
        nName = seq.replaceName(sName, namedict, "rename")
        lens[nName] = len(sSeq)
        names[nName] = sName.replace("piece", "tig")
        #sys.stderr.write("Updating name to be mapping %s to %s of len %s\n"%(nName, sName.replace("piece", "tig"), len(sSeq)))
     else:
        sys.stderr.write("Error: unexpected input %s\n"%(line))
        sys.exit(1)

# save the offets of individual tigs and also update the header to include our new references
header = bamfile.header.to_dict()

for clist in scfmap:
   offset = 0
   for piece in scfmap[clist]:
      numn = re.match(r"\[N(\d+)N]", piece)
      if numn:
         offset += int(numn[1])
      elif piece in lens:
         offsets[names[piece]] = offset
         #sys.stderr.write("The offset for %s is %s\n"%(names[piece], offset))
         offset += lens[piece]
         tigtohap[names[piece]] = clist
   #sys.stderr.write("Saving len of %s for sequence %s\n"%(offset, clist))     
   header['SQ'].append({'SN': clist, 'LN': offset})

# drop the old reference names
header['SQ'] = [sq_entry for sq_entry in header['SQ'] if sq_entry['SN'] in scfmap]

# add read tags to note which reads are LA, LA gapfill, or UL
if 'RG' not in header:
    header['RG'] = []
header['RG'].append({'ID' : "LA"})
header['RG'].append({'ID' : "UL-gap"})
header['RG'].append({'ID' : "UL"})

# now loop through and update, keeping only the mapped reads
output_bam = pysam.AlignmentFile('-', 'wb', header=header)
for read in bamfile:
    if not read.is_unmapped:
       reference_name = bamfile.get_reference_name(read.reference_id)
       if reference_name not in offsets or tigtohap[reference_name] not in output_bam.references or read.query_name not in readtorg:
          sys.stderr.write("Error: I didn't find how to translate '%s'\n"%(reference_name))
          sys.exit(1)
     
       # make a copy of the existing read in the new bam, replacing the string reference
       new_read = pysam.AlignedSegment.fromstring(read.to_string().replace(f'\t{reference_name}\t', f'\t{tigtohap[reference_name]}\t'), output_bam.header)

       new_read.set_tag('RG', readtorg[read.query_name])
       # update the name to strip layout prefix
       new_read.query_name = next((new_read.query_name[len(p):] for p in prefixes if new_read.query_name.startswith(p)), new_read.query_name)
       # update the coordinate
       new_read.reference_start = read.reference_start + offsets[reference_name]
       output_bam.write(new_read)

bamfile.close()
output_bam.close()
