#!/usr/bin/env python3

import sys
import re
import os

def print_results(names):
   for i in range(len(names)):
      for j in range(i+1, len(names)):
         i_split = names[i].split()
         j_split = names[j].split()
         if (i_split[0] != j_split[0]):
            print("Error: name read %s and %s don't match"%(i_split, j_split), file=sys.stderr)
            sys.exit(1)
         if (i_split[1] < j_split[1]):
            print (f"{i_split[0]}\t{i_split[1]}\t{j_split[1]}\t1\t{i_split[2]}\t{j_split[2]}")
         elif (i_split[1] > j_split[1]):
            print (f"{i_split[0]}\t{j_split[1]}\t{i_split[1]}\t1\t{j_split[2]}\t{i_split[2]}")
         # if they are equal we don't want to print anything, not informative 

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
        input_stream = open(input_filename, 'r')

name = ""
names = [ ]

for line in input_stream:
   line=line.split()
   if len(line) < 4:
      continue
   if name == "":
      name = line[0]
   if name != line[0]:
      print_results(names)
      name = line[0]
      names = [ ]
   if name in seen:
      print("Warning: read %s already seen but encountered it again, please confirm your bam file is sorted by read."%(name), file=sys.stderr)
      out_of_order += 1
   names.append("%s\t%s\t%s"%(line[0], line[2], line[3]))

print_results(names)
