import sys

layout_file = sys.argv[1]
aln_file = sys.argv[2]

def process_alns(aln_file):
        aln_dict = {}
        alns = open(aln_file)
        for line in alns:
                line_split = line.strip().split('\t')
                qname, start, end, orient, rname = line_split[0], int(line_split[2]), int(line_split[3]), line_split[4], line_split[5]
                k = qname + '_to_' + rname
                if k in alns:
                        alns[k] += {start, end, orient}}
                else:
                        alns[k] = [{start, end, orient}]]
        return alns_dict


aln_dict = process_alns(aln_file)

for line in layout_file:
        if line.startswith("tig"):
                current_tig = line.split("\t")[1]
                num_read_seen = dict()
        elif line.startswith("len"):
                length = int(line.split("\t")[1])
        else:
                line_split = line.strip().split('\t')
                read_name, start, end, is_ont = line_split[0], int(line_split[1]), int(line_split[2]), line_split[3]
                if read_name in num_read_seen:
                        num_read_seen[read_name] += 1
                else:
                        num_read_seen[read_name] = 1
                if is_ont and end > length:
                        k = read_name + '_to_' + current_tig
                        entry = aln_dict[k][num_read_seen[read_name]]
                        new_line = line.strip() + '\t' + entry['start'] + '\t' + entry[end]
                        print(new_line)

