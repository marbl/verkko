import sys

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
            elif line == "end":
                current_tig = None
            else:
                read, start, end = line.split("\t")
                layout_data[-1]["reads"].append({"read": read, "start": int(start), "end": int(end), "flag": flag})

    return layout_data


def merge_layout_files(file1, flag1, file2, flag2):
    layout1 = read_layout_file(file1, flag1)
    layout2 = read_layout_file(file2, flag2)

    merged_layout = []

    # Merge layout data for tigs present in both files
    for layout1_entry in layout1:
        tig1 = layout1_entry["tig"]
        layout2_entry_found = False

        for layout2_entry in layout2:
            if layout2_entry["tig"] == tig1:
                layout2_entry_fixed_reads = []
                layout2_entry_fixed_num_reads = 0
                for read in layout2_entry["reads"]:
                    if read["start"] >= layout1_entry["len"] and read["end"] >= layout1_entry["len"]:
                        continue
                    # if read["end"] >= layout1_entry["len"]:
                    #     read["end"] = layout1_entry["len"]
                    layout2_entry_fixed_reads.append(read)
                    layout2_entry_fixed_num_reads += 1
                merged_entry = {
                    "tig": tig1,
                    "len": layout1_entry["len"],
                    "num_reads": layout1_entry["num_reads"] + layout2_entry_fixed_num_reads,
                    "reads": layout1_entry["reads"] + layout2_entry_fixed_reads,
                }
                merged_entry["reads"].sort(key=lambda x: min(x['start'], x['end']))
                layout2.remove(layout2_entry)
                layout2_entry_found = True
                break

        if not layout2_entry_found:
            merged_entry = {
                "tig": tig1,
                "len": layout1_entry["len"],
                "num_reads": layout1_entry["num_reads"],
                "reads": layout1_entry["reads"],
            }

        merged_layout.append(merged_entry)

    # Add remaining layout data from layout2
    merged_layout.extend(layout2)

    return merged_layout


file1_path = sys.argv[1]
file2_path = sys.argv[2]

merged_layout = merge_layout_files(file1_path, 0, file2_path, 1)

# Print the merged layout
for entry in merged_layout:
    print("tig", entry["tig"], sep='\t')
    print("len", entry["len"], sep='\t')
    print("rds", entry["num_reads"], sep='\t')
    
    for read_entry in entry["reads"]:
        print(read_entry["read"], read_entry["start"], read_entry["end"], read_entry['flag'], sep='\t')
    
    print("end")
