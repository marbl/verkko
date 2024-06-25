import sys

def utig1_to_utig3(node, pos, ref, alt):
    return

def utig3_to_utig4(node, pos, ref, alt):
    return

def utig4_to_rukki_path(node, pos, ref, alt):
    return

def rukki_path_to_output_hap(node, pos, ref, alt):
    return

def utig1_to_output_hap(node, pos, ref, alt):
    return rukki_path_to_output_hap(utig4_to_rukki_path(utig3_to_utig4(utig1_to_utig3(node, pos, ref, alt))))


in_vcf = open(sys.argv[1])
verkko_folder = sys.argv[2]
map_3 = open(verkko_folder + '/5-untip/unitig-mapping-3.txt')
map_4 = open(verkko_folder + '/5-untip/unitig-mapping-4.txt')
unitig_paths_gaf = open(verkko_folder + '/6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf')
unitig_layout = open(verkko_folder + '/6-layoutContigs/unitig-popped.layout.scfmap')

for line in in_vcf:
    line_split = line.strip().split('\t')
    og_node, og_pos, og_ref, og_alt = line_split[0], line_split[1], line_split[3], line_split[4]
    final_node, final_pos, final_ref, final_alt = utig1_to_output_hap(og_node, og_pos, og_ref, og_alt)
    new_line = '\t'.join([final_node, final_pos] + [line_split[2]] + [final_ref, final_alt] + line_split[5:])
    print(new_line)
