#!/usr/bin/env python

import sys

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
                read, start, end = line.split("\t")
                layout_data[-1]["reads"].append({"read": read, "start": int(start), "end": int(end), "flag": flag})

    return layout_data


def merge_layout_files(file1, flag1, file2, flag2):
    layout1 = read_layout_file(file1, flag1)
    layout2 = read_list(file2)

    merged_layout = []

    # Merge layout data for tigs present in both files
    for layout1_entry in layout1:
        keepTig = False
        tig1 = layout1_entry["tig"]

        for read in layout1_entry["reads"]:
           if read['read'] in layout2:
              read['flag'] = flag2
           else:
              keepTig = True

        if keepTig == True: merged_layout.append(layout1_entry)

    return merged_layout


file1_path = sys.argv[1]
file2_path = sys.argv[2]

merged_layout = merge_layout_files(file1_path, 0, file2_path, 1)

# Print the merged layout
for entry in merged_layout:
    print("tig", entry["tig"], sep='\t')
    print("len", entry["len"], sep='\t')
    print("trm", entry["trm"], sep='\t')
    print("rds", entry["num_reads"], sep='\t')

    for read_entry in entry["reads"]:
        print(read_entry["read"], read_entry["start"], read_entry["end"], read_entry['flag'], sep='\t')

    print("end")
