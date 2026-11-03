#!/usr/bin/env python
import pysam
import sys
import fasta_util as seq
import re

# Rename alignments based on the input scfmap and tig name mapping files.
# CRAM inputs are opened sequentially so package count does not translate into
# open-file pressure.  With no --alignments argument, retain the legacy BAM
# stdin interface.
# in the case of scaffolds, offset the start coordinate by the start of the sequence in the scaffold
#
prefixes = [ 'hifi_', 'ont_' ]
layout   = sys.argv[1]
gap_info = sys.argv[2]
scfmap,_ = seq.readScfMap(sys.argv[3])
namedict = seq.readNameMap(sys.argv[4])

if '--alignments' in sys.argv:
   separator = sys.argv.index('--alignments')
   consensus_files = sys.argv[5:separator]
   alignment_files = sys.argv[separator + 1:]
   if not alignment_files:
      sys.stderr.write("Error: --alignments requires at least one CRAM file\n")
      sys.exit(1)
else:
   consensus_files = sys.argv[5:]
   alignment_files = []

lens     = dict()
offsets  = dict()
names    = dict()
tigtohap = dict()
readtorg = dict()

with open(layout) as f:
   for l in f:
      parts = l.strip().split('\t')
      if len(parts) > 2:
         readtorg[parts[0]] = "LA" if (int(parts[-1]) == 0 or len(parts) == 4) else "UL"
with open(gap_info) as f:
   for l in f:
      parts = l.strip().split('\t')
      if parts[0] in readtorg and readtorg[parts[0]] == "LA":
          readtorg[parts[0]] = "UL-gap"

for filename in consensus_files:
  sys.stderr.write("Starting file %s\n"%(filename))
  inf  = seq.openInput(filename)

  line = inf.readline()

  while (line != ""):
     if   (line[0] == ">"):
        line, sName, sSeq, sQlt = seq.readFastA(inf, line)
        nName = seq.replaceName(sName, namedict, "rename")
        if len(sSeq) > 0:
           lens[nName] = len(sSeq)
           names[nName] = sName.replace("piece", "tig")
           #sys.stderr.write("Updating name to be mapping %s to %s of len %s\n"%(nName, sName.replace("piece", "tig"), len(sSeq)))
     else:
        sys.stderr.write("Error: unexpected input %s\n"%(line))
        sys.exit(1)

# Save the offsets of individual tigs and update the header with the new
# references.  Old package SQ entries are discarded below; read groups are
# unioned across package headers so caller tier metadata is retained.
if alignment_files:
   input_read_groups = []
   seen_read_groups = set()
   header = None
   # Inspect headers sequentially too.  Package CRAMs normally have identical
   # tier groups, but unioning them avoids losing a group present only in a
   # later package without increasing the open-file count.
   for alignment_file in alignment_files:
      with pysam.AlignmentFile(alignment_file, "rc") as header_file:
         input_header = header_file.header.to_dict()
         if header is None:
            header = input_header
         for group in input_header.get('RG', []):
            group_id = group.get('ID')
            if group_id is not None and group_id not in seen_read_groups:
               input_read_groups.append(group)
               seen_read_groups.add(group_id)
   bamfile = None
else:
   bamfile = pysam.AlignmentFile("-", "rb")
   header = bamfile.header.to_dict()
   input_read_groups = header.get('RG', [])

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
   if (offset > 0):
      #sys.stderr.write("Saving len of %s for sequence %s\n"%(offset, clist))
      header['SQ'].append({'SN': clist, 'LN': offset})

# drop the old reference names
header['SQ'] = [sq_entry for sq_entry in header['SQ'] if sq_entry['SN'] in scfmap]

# Preserve caller-supplied groups such as Loimi's LA/LAO/UL tiers.  Gap-filling
# reads are the only records reclassified in this mode.  Legacy BAMs without
# groups retain the original layout-derived LA/UL/UL-gap assignment.
preserve_read_groups = bool(input_read_groups)
if preserve_read_groups:
   header['RG'] = input_read_groups
   if not any(group.get('ID') == 'UL-gap' for group in header['RG']):
      header['RG'].append({'ID': 'UL-gap'})
else:
   header['RG'] = [
      {'ID': "LA"},
      {'ID': "UL-gap"},
      {'ID': "UL"},
   ]

# Now loop through and update, keeping only the mapped reads.  The rewritten
# BAM stream is transient; the workflow immediately sorts and encodes it as
# the final embedded-reference CRAM.
output_bam = pysam.AlignmentFile('-', 'wb', header=header)

def rename_alignments(input_file):
   for read in input_file:
      if not read.is_unmapped:
         reference_name = input_file.get_reference_name(read.reference_id)
         if reference_name not in offsets or tigtohap[reference_name] not in output_bam.references or read.query_name not in readtorg:
            sys.stderr.write("Error: I didn't find how to translate '%s' hap is '%s' and query name '%s'\n"%(reference_name, tigtohap[reference_name], read.query_name))
            sys.exit(1)

         # Make a copy of the existing read using the new reference header.
         new_read = pysam.AlignedSegment.fromstring(read.to_string().replace(f'\t{reference_name}\t', f'\t{tigtohap[reference_name]}\t'), output_bam.header)

         if readtorg[read.query_name] == 'UL-gap':
            new_read.set_tag('RG', 'UL-gap')
         elif not preserve_read_groups:
            new_read.set_tag('RG', readtorg[read.query_name])
         # Update the name to strip layout prefix.
         new_read.query_name = next((new_read.query_name[len(p):] for p in prefixes if new_read.query_name.startswith(p)), new_read.query_name)
         # Lift the piece-relative coordinate into the scaffold.
         new_read.reference_start = read.reference_start + offsets[reference_name]
         output_bam.write(new_read)

if alignment_files:
   for alignment_file in alignment_files:
      sys.stderr.write("Reading alignments from %s\n" % alignment_file)
      with pysam.AlignmentFile(alignment_file, "rc") as input_file:
         rename_alignments(input_file)
else:
   rename_alignments(bamfile)
   bamfile.close()

output_bam.close()
