#!/usr/bin/env python
import pysam
import sys

lens = {}
cur_name = ""
for line in open ("unitigs.fasta"):
    if line[0] == ">":
        cur_name = line[1:].strip()
    else: 
        lens[cur_name] = int(len(line.strip()))
bamfile = pysam.AlignmentFile("-", "rb")

total_reads = 0
unique_pairs = 0
valid_pairs = 0
all_pairs = 0
used_pairs = 0
cur_name = ""
reads = []
prev_read = None
prev_name = ""
unpaired_reads = 0

both_short_total = 0
one_short_total = 0
both_middle_total = 0
one_middle_total = 0
basic_filtered = 0

XA_self_mapped = 0
output_bam = pysam.AlignmentFile('-', 'wb', template=bamfile)
for cur_read in bamfile:
    total_reads += 1
# awk '$7 != "*" && $7 != "="' | awk '$NF !~ "SA" && $NF !~ ($7 "," )' | awk '$5 != 0 || $NF ~ "XA"'
    if cur_read.is_unmapped or cur_read.reference_name == cur_read.next_reference_name or (cur_read.has_tag("XA") and cur_read.get_tag("XA").find(cur_read.next_reference_name+",") != -1) or (not cur_read.has_tag("XA") and cur_read.mapping_quality == 0):
        basic_filtered += 1
        continue
    cur_name = cur_read.query_name
    if cur_name == prev_name:
        #TODO: poreC is not compatible with this now
        #Last read is always missing but who cares? do not want to make it function since it is time-critical part
        reads = (prev_read, cur_read)
#                  if read.is_paired:
        all_pairs += 1
        unique = False
        if prev_read.mapping_quality > 0 and cur_read.mapping_quality > 0:
            #TODO: special storage
            unique_pairs += 1
            unique = True
            
        names = [[prev_read.reference_name], [cur_read.reference_name]]
        coords = [[prev_read.reference_start], [cur_read.reference_start]]
#                    names = read.reference_name
        valid = True
        i = 0
        for read in reads:
            if read.has_tag("XA"):
                for xa in read.get_tag("XA")[:-1].split(";"):
                    xa_arr = xa.split(",")
                    names[i].append(xa_arr[0])
                    coords[i].append(int(xa_arr[1][1:]))
            i += 1
        self_mapped = False
        for name0 in names[0]:
            for name1 in names[1]:
                if name0 == name1:
                    self_mapped = True
                    XA_self_mapped = True
        all_short = [True, True]
        all_middle = [True, True]
        if not self_mapped:
            for ind in range (0, 2):
                #print (names[ind])
                #print ( len(names[ind]) )
                for i in range (0, len(names[ind]) ):
                    node_f_len = lens[names[ind][i]]
                    if node_f_len > 20000: 
                        all_short[ind] = False
                    if coords[ind][i] < 10000000 and node_f_len - coords[ind][i] < 10000000:
                        all_middle[ind] = False
            if all_short[0] and all_short[1]:
                both_short_total += 1
            if all_short[0] or all_short[1]:
                one_short_total += 1
            if all_middle[0] and all_middle[1]:
                both_middle_total += 1
            if all_middle[0] or all_middle[1]:
                one_middle_total += 1
        if unique or ((not (all_short[0] or all_short[1])) and (not (all_middle[0] or all_middle[1])) and not self_mapped):
            output_bam.write(prev_read)
            output_bam.write(cur_read)

    else:
        unpaired_reads += 1

    prev_read = cur_read
    prev_name = cur_name
bamfile.close()
output_bam.close()

sys.stderr.write (f"total reads {total_reads} unique pairs {unique_pairs} paired {all_pairs}\n")
sys.stderr.write (f"Not paired {unpaired_reads - all_pairs} basic filtered {basic_filtered} \n")
sys.stderr.write (f"For pairs addfiltered: Both short {both_short_total} one short {one_short_total} both middle {both_middle_total} one middle {one_middle_total} XA_self_mapped {XA_self_mapped}\n")