#!/usr/bin/env python

import sys
import fasta_util as seq
import re

def read_list(file_path):
   read_list = set()
   with open(file_path, "r") as file:
       for line in file:
           line = line.strip().split("\t")
           read_list.add(line[0])
   return read_list

def read_layout_file(file_path, flag):
    layout_data = []
    current_tig = None

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if line.startswith("tig"):
                current_tig = line.split("\t")[1]
            elif line.startswith("len"):
                length = int(line.split("\t")[1])
                layout_data.append({"tig": current_tig, "len": length, "reads": []})
            elif line.startswith("rds"):
                num_reads = int(line.split("\t")[1])
                layout_data[-1]["num_reads"] = num_reads
            elif line.startswith("trm"):
                layout_data[-1]["trm"] = int(line.split("\t")[1])
            elif line == "end":
                current_tig = None
            else:
                read, start, end, copy = line.split("\t")
                layout_data[-1]["reads"].append({"read": read, "start": int(start), "end": int(end), "copy" : int(copy), "flag": flag})

    return layout_data


def merge_layout_files(file1, flag1, file2, flag2):
    layout1 = read_layout_file(file1, flag1)
    layout2 = read_list(file2)

    merged_layout = []
    dropped = set()

    # Merge layout data for tigs present in both files
    for layout1_entry in layout1:
        keepTig = False
        tig1 = layout1_entry["tig"]

        for read in layout1_entry["reads"]:
           if read['read'] in layout2:
              read['flag'] = flag2
           else:
              keepTig = True

        if keepTig == True:
            merged_layout.append(layout1_entry)
        else:
            dropped.add(layout1_entry["tig"])

    return merged_layout, dropped


file1_path      = sys.argv[1]
file2_path      = sys.argv[2]
output_prefix   = sys.argv[3]
scfmap, scfname = seq.readScfMap(f"{sys.argv[1]}.scfmap")
tig_layout_file = open(f"{output_prefix}", mode="w")
scf_layout_file = open(f"{output_prefix}.scfmap", mode="w")
nul_layout_file = open(f"{output_prefix}.scfmap.dropped", mode="a")

merged_layout, dropped = merge_layout_files(file1_path, 0, file2_path, 1)

# Print the merged layout
for entry in merged_layout:
    print("tig", entry["tig"], sep='\t', file=tig_layout_file)
    print("len", entry["len"], sep='\t', file=tig_layout_file)
    print("trm", entry["trm"], sep='\t', file=tig_layout_file)
    print("rds", entry["num_reads"], sep='\t', file=tig_layout_file)

    for read_entry in entry["reads"]:
        print(read_entry["read"], read_entry["start"], read_entry["end"], read_entry['copy'], read_entry['flag'], sep='\t', file=tig_layout_file)

    print("end", file=tig_layout_file)

# print the merged scfmap
for clist in scfmap:
   keep = False
   for piece in scfmap[clist]:
      numn = re.match(r"\[N(\d+)N]", piece)
      if numn or piece not in dropped:
         keep= True

   # now that we've checked that some pieces are saved, output it
   if keep == False:
      print(f"{clist} has no reads assigned and is not output.", file=nul_layout_file)
      continue

   previous = "NONE"
   new_pieces = []
   print(f"path {clist} {scfname[clist]}", file=scf_layout_file)
   for piece in scfmap[clist]:
      is_gap = re.match(r"\[N\d+N\]", piece)
      if is_gap:
         if previous == "NONE":
            continue
         if previous == "piece":
            new_pieces.append(f"{piece}")
            previous = "gap"
         else:
            #merging consecutive gaps
            last = new_pieces.pop()
            previous = "gap"
            last_int = int(last[2:-2])
            cur_int = int(piece[2:-2])
            new_pieces.append("[N%dN]"%(last_int + cur_int))
      else:
         if piece in dropped: continue
         new_pieces.append(f"{piece}")
         previous = "piece"
   print("\n".join(new_pieces), file=scf_layout_file)
   print("end", file=scf_layout_file)
