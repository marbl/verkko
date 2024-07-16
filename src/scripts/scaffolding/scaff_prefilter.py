#!/usr/bin/env python
import pysam
import sys


total_reads = 0
unique_pairs = 0
valid_pairs = 0
all_pairs = 0
used_pairs = 0
unpaired_reads = 0

both_short_total = 0
one_short_total = 0
both_middle_total = 0
one_middle_total = 0
basic_filtered = 0

XA_self_mapped = 0

bamfile = pysam.AlignmentFile(sys.argv[1], "rb")
output_bam = pysam.AlignmentFile(sys.argv[2], 'wb', template=bamfile)
interesting_nodes = set()
for line in open (sys.argv[3]):
    interesting_nodes.add(line.strip())

cur_name = ""
reads = []
prev_read = None
prev_name = ""

for cur_read in bamfile:
    total_reads += 1
    cur_name = cur_read.query_name
    if cur_name == prev_name:
        #TODO: poreC is not compatible with this now
        reads = (prev_read, cur_read)
        all_pairs += 1
        if prev_read.mapping_quality > 0 and cur_read.mapping_quality > 0:
            unique_pairs += 1
        names = [[prev_read.reference_name], [cur_read.reference_name]]
        coords = [[prev_read.reference_start], [cur_read.reference_start]]
        quals = [prev_read.mapping_quality, cur_read.mapping_quality]
        valid = True
        i = 0
        for read in reads:
            if read.has_tag("XA"):
                nm = int(read.get_tag("NM"))
                if read.has_tag("XS") and read.get_tag("XS") == read.get_tag("AS"):  
                    for xa in read.get_tag("XA")[:-1].split(";"):
                        xa_arr = xa.split(",")
                        #do not want to do cigar comparsion
                        if int(xa_arr[3]) == nm:
                            names[i].append(xa_arr[0])
            i += 1
        filtered_names = [[], []]
        filtered_coords = [[], []]
        for i in range(0, 2):
            lname = len(names[i])
            for j in  range (0, lname):                
                if names[i][j] in interesting_nodes:
                    filtered_names[i].append(names[i][j])

        lname0 = len(filtered_names[0])
        lname1 = len(filtered_names[1])
        #self.logger.debug (f" {valid} {lname0} {lname1}")
        if valid and lname0 > 0 and lname1 > 0:
            valid_pairs += 1
            output_bam.write(reads[0])
            output_bam.write(reads[1])
        else:
            basic_filtered += 1
    prev_read = cur_read
    prev_name = cur_name
bamfile.close()
output_bam.close()

sys.stderr.write (f"total read_pairs {total_reads/2} filtering with list {basic_filtered}\n")
