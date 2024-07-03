#!/usr/bin/env python3

import sys
import re
import os

def count_matches(cigar_string):
    # Find all the occurrences of the pattern in the CIGAR string
    
    matches = re.findall(r'(\d+)M', cigar_string)
    # Convert the list of strings to integers and sum them up
    total_matches = sum(map(int, matches))
    #print (f"{cigar_string} - {total_matches}", file=sys.stderr)
    return total_matches

def print_results(names):
   nlen = len(names)
   for i in range(0, nlen - 1):
      for j in range(i+1, nlen):
         aln_split = [names[i], names[j]]

         if (aln_split[0][0] != aln_split[1][0]):
            print("Error: name read %s and %s don't match"%(aln_split[0], aln_split[1]), file=sys.stderr)
            sys.exit(1)

         if (aln_split[0][2] > aln_split[1][2]):
            aln_split.reverse()
         valid = True
         for aln in aln_split:
            if aln[2] == "*":
               valid = False
         if not valid:
            continue
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
            # if we do not have alternative alignments and mapping quality is 0 - skip, too bad or too much alternatives
            if not ("XA" in tags[i]) and aln_split[i][4] == "0":
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
            primary_matches = count_matches(sams[i])
            if "XA" in tags[i]:
               for alt in tags[i].get("XA").split(";")[:-1]:
                  alt_split = alt.split(",")
                  #TODO: should we really forbid it?
                  if alt_split[0] in matched_ids[1 - i]:
                     valid = False
                     break
                  #some of alts can be same as primary, some can be worse
                  if count_matches(alt_split[2]) == primary_matches and int(alt_split[3]) <= int(tags[i]["NM"]):
                     matched_ids[i].append(alt_split[0])
                     matched_coords[i].append(alt_split[1][1:])
            if len(matched_ids[i]) >= 5:
               valid = False
         if not valid:
            continue
         #if we have alignment to the same contig, then it's likely the answer and no multimappers needed
         #should be faster than sets comparisons for small arrays
         for first_id in matched_ids[0]:
            for second_id in matched_ids[1]:
               if first_id == second_id:
                  valid = False
         if not valid:
            continue
         #print (f'{aln_split[0][0]}\t{",".join(matched_ids[0])}\t{",".join(matched_ids[1])}\t1\t{",".join(matched_coords[0])}\t{",".join(matched_coords[1])}')
         print (f'X\t{",".join(matched_ids[0])}\t{",".join(matched_ids[1])}\t1\t{",".join(matched_coords[0])}\t{",".join(matched_coords[1])}')

         
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
