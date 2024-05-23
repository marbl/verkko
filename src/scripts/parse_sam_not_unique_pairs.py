#!/usr/bin/env python3

import sys
import re
import os

def print_results(names):
   for i in range(len(names)):
      for j in range(i+1, len(names)):
         aln_split = [names[i].split(), names[j].split()]

         if (aln_split[0][0] != aln_split[1][0]):
            print("Error: name read %s and %s don't match"%(aln_split[0], aln_split[1]), file=sys.stderr)
            sys.exit(1)
         if (aln_split[0][1] > aln_split[1][1]):
            aln_split.reverse()
         valid = True
         for aln in aln_split:
            if aln[2] == "*":
               valid = False
         tags = [{},{}]
         for i in range(0, 2):
            for tag in aln_split[i][11:]:
               tag_split = tag.split(":")
               if len(tag_split) < 3:
                  print ("invalid tag %s"%(tag), file=sys.stderr)
                  sys.exit(1)
               tags[i][tag_split[0]] = tag_split[2]
         sams = [aln_split[0][5], aln_split[1][5]]
         for i in range(0, 2):
            # if we do not have alternative alignments and mapping quality is 0 - skip
            if not ("XA" in tags[i]) and aln_split[4] == "0":
               valid = False
            #need sam and edit distance to compare reported second bests. Not using XS because of multiple suboptimals
            if not ("NM" in tags[i]):
               valid = False
         if not valid:
            continue
         matched_ids = [[aln_split[0][2]], [aln_split[1][2]]]
         matched_coords = [[aln_split[0][3]], [aln_split[1][3]]]
         #XA:Z:utig4-2161,-5673767,151M,3;utig4-183,+3013486,151M,4;utig4-184,-3355570,151M,5;utig4-184,+4936592,151M,5;
         for i in range (0, 2):
            if "XA" in tags[i]:
               for alt in tags[i].get("XA").split(";")[:-1]:
                  alt_split = alt.split(",")
                  #some of alts can be same as primary, some can be worse
                  if alt_split[2] == sams[i] and int(alt_split[3]) <= int(tags[i]["NM"]):
                     matched_ids[i].append(alt_split[0])
                     matched_coords[i].append(alt_split[1][1:])
            if len(matched_ids[i]) > 5:
               valid = False
         if not valid:
            continue
         #if we have alignment to the same contig, then it's likely the answer and no multimappers needed
         for first_id in matched_ids[0]:
            for second_id in matched_ids[1]:
               if first_id == second_id:
                  valid = False
         if not valid:
            continue
         len0 = len(matched_ids[0])
         len1 = len(matched_ids[1])
         inv_weight = len0 * len1
         for i in range(0, len0):
            for j in range(0, len1):
               ids = [matched_ids[0][i], matched_ids[1][j]]
               coords = [matched_coords[0][i], matched_coords[1][j]]
               if ids[0] > ids[1]:
                  ids.reverse()
                  coords.reverse()
               print (f"{aln_split[0][0]}\t{ids[0]}\t{ids[1]}\t{inv_weight}\t{coords[0]}\t{coords[1]}")

if not sys.stdin.isatty():
    input_stream = sys.stdin

# otherwise, read the given filename                                            
else:
    try:
        input_filename = sys.argv[1]
    except IndexError:
        message = 'need filename as first argument if stdin is not full'
        raise IndexError(message)
    else:
        input_stream = open(input_filename, 'rU')

name = ""
names = [ ]
seen = {}
out_of_order = 0

for line in input_stream:
   line=line.split()
   if name == "":
      name = line[0]
   if name != line[0]:
      seen[name] = 1
      print_results(names)
      name = line[0]
      names = [ ]
   if name in seen:
      print("Warning: read %s already seen but encountered it again, please confirm your bam file is sorted by read."%(name), file=sys.stderr)
      out_of_order += 1
   names.append(line)

if out_of_order > 1000:
   print("Error: encountered too many unsorted reads (%d), exiting. Please confirm the input bam is sorted by read."%(out_of_order), file=sys.stderr)
   sys.exit(1)

print_results(names)
