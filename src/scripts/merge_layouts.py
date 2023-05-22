import re
import sys

hifi_layout = open(sys.argv[1])
ont_layout = open(sys.argv[2])

lines = dict()
for hifi_line in hifi_layout:
    if hifi_line[:3] == 'tig':
        hifi_tig_line = hifi_line.strip()
    elif hifi_line[:3] == 'len':
        hifi_len_line = hifi_line.strip()
    elif hifi_line[:3] == 'rds':
        hifi_num_reads = int(hifi_line.strip().split('\t')[1])
    elif re.match('^[a-zA-Z0-9-_\/]+\t\d+\t\d+$', hifi_line):
        hifi_line_split = hifi_line.split('\t')
        hifi_pos = min(int(hifi_line_split[1]), int(hifi_line_split[2]))
        lines[hifi_line.strip() + '\t0'] = hifi_pos
    elif hifi_line[:3] == 'end':
        ont_num_reads = -1
        keep_lines = False
        printed_header_and_hifi = False
        for ont_line in ont_layout:
            if ont_line.strip() == hifi_tig_line:
                keep_lines = True
            if keep_lines and ont_line[:3] == 'rds':
                ont_num_reads = int(ont_line.strip().split('\t')[1])
            if keep_lines and ont_num_reads >= 0 and not printed_header_and_hifi:
                print(hifi_tig_line)
                print(hifi_len_line)
                print('rds\t' + str(hifi_num_reads + ont_num_reads))
                printed_header_and_hifi = True
            if keep_lines and re.match('^[a-zA-Z0-9-_\/]+\t\d+\t\d+$', ont_line):
                ont_line_split = ont_line.split('\t')
                ont_pos = min(int(ont_line_split[1]), int(ont_line_split[2]))
                lines[ont_line.strip() + '\t1'] = ont_pos
            if ont_line[:3] == 'end':
                sorted_lines = dict(sorted(lines.items(), key=lambda x:x[1]))
                for k,v in sorted_lines.items():
                    print(k)
                print('end')
                hifi_num_reads = -1
                lines = dict()
                break


