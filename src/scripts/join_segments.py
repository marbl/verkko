#!/usr/bin/env python
from __future__ import division
from __future__ import print_function

import argparse
import re
import sys


def normalize_segment_name(name):
    """Match the historical Verkko convention of dropping the tig zero prefix."""
    return re.sub(r"tig0+", "", name)


def my_segment(sequence, orientation):
    if orientation == "-":
        return sequence.reverse_complement()
    if orientation == "+":
        return sequence
    raise ValueError("Unknown segment orientation %r" % orientation)


def swap(orientation):
    if orientation == "-":
        return "+"
    if orientation == "+":
        return "-"
    raise ValueError("Unknown segment orientation %r" % orientation)


def read_graph(gfa_fn):
    """Read segment sequences and bidirectional link overlaps from a GFA."""
    from Bio.Seq import Seq

    segments = {}
    cigars = {}
    with open(gfa_fn) as gfa:
        for line in gfa:
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue

            if fields[0] == "S":
                if len(fields) < 3:
                    raise ValueError("GFA segment has fewer than three fields: %s" % line.rstrip())
                if fields[2] != "*":
                    segments[normalize_segment_name(fields[1])] = Seq(fields[2])
                continue

            if fields[0] != "L":
                continue
            if len(fields) < 6:
                raise ValueError("GFA link has fewer than six fields: %s" % line.rstrip())

            first = normalize_segment_name(fields[1]) + fields[2]
            second = normalize_segment_name(fields[3]) + fields[4]
            cigars[first + "|" + second] = fields[5]
            reverse = normalize_segment_name(fields[3]) + swap(fields[4])
            reverse_next = normalize_segment_name(fields[1]) + swap(fields[2])
            cigars[reverse + "|" + reverse_next] = fields[5]
    return segments, cigars


def get_cigar(cigars, first, second):
    return cigars.get(first + "|" + second)


def trim_len(cigar, second=True):
    """Return the sequence consumed on the requested side of a link CIGAR."""
    pos = sum(int(match.group(1)) for match in re.finditer(r"(\d+)M", cigar))
    indel = r"(\d+)I" if second else r"(\d+)D"
    pos += sum(int(match.group(1)) for match in re.finditer(indel, cigar))
    return pos


def parse_gaf_path(path):
    """Parse a contiguous Rukki GAF-style path into oriented segments.

    Segment names extend up to the next orientation marker or gap, so names such
    as ``utig4-443`` and ``node.1`` are preserved. Gapped paths are intentionally
    rejected because gaps are resolved before this script is run.
    """
    if not path:
        raise ValueError("Path is empty")

    steps = []
    pos = 0
    while pos < len(path):
        marker = path[pos]
        if marker in "><":
            end = pos + 1
            while end < len(path) and path[end] not in "><[":
                end += 1
            name = path[pos + 1:end]
            if not name:
                raise ValueError("Missing segment name at offset %d in %r" % (pos, path))
            orientation = "+" if marker == ">" else "-"
            steps.append(("segment", normalize_segment_name(name), orientation))
            pos = end
            continue

        if marker == "[":
            raise ValueError("Gapped paths are not supported: %r" % path)

        raise ValueError("Unexpected character at offset %d in %r" % (pos, path))

    return steps


def make_sequence(cigars, contig_dict, steps, trim_flanks_to):
    """Spell a parsed path, trimming the graph overlap at each link."""
    first_desc = steps[0]
    first_name, first_orientation = first_desc[1], first_desc[2]
    first_sequence = my_segment(contig_dict[first_name], first_orientation)

    if trim_flanks_to < 0:
        final = first_sequence
    else:
        final = first_sequence[-trim_flanks_to:]
        print("Here taking", min(trim_flanks_to, len(first_sequence)), file=sys.stderr)

    previous_desc = first_desc
    for step in steps[1:]:
        name, orientation = step[1], step[2]
        sequence = my_segment(contig_dict[name], orientation)
        previous = previous_desc[1] + previous_desc[2]
        current = name + orientation
        cigar = get_cigar(cigars, previous, current)
        if cigar is None:
            raise ValueError("CIGAR not known between %s and %s" % (previous, current))
        final += sequence[trim_len(cigar):]
        previous_desc = step

    last_sequence = my_segment(contig_dict[steps[-1][1]], steps[-1][2])
    if trim_flanks_to > 0 and len(last_sequence) > trim_flanks_to:
        final = final[:-(len(last_sequence) - trim_flanks_to)]

    return final, False


def join_paths(gfa_path, paths_path, trim_flanks_to):
    """Yield spelled FASTA records in path order without buffering outputs."""
    from Bio.SeqRecord import SeqRecord

    print("Reading sequences and links from", gfa_path, file=sys.stderr)
    contig_dict, cigars = read_graph(gfa_path)

    with open(paths_path) as paths:
        for line_number, line in enumerate(paths, 1):
            fields = line.split()
            if not fields or fields[0] == "name":
                continue
            if len(fields) < 2:
                raise ValueError("Path line %d has fewer than two fields" % line_number)

            name, path = fields[0], fields[1]
            try:
                steps = parse_gaf_path(path)
                sequence, unused_circular = make_sequence(
                    cigars, contig_dict, steps, trim_flanks_to
                )
            except (KeyError, ValueError) as error:
                raise ValueError(
                    "Failed to join path %s on line %d: %s" % (name, line_number, error)
                ) from error
            yield SeqRecord(sequence, id=name, description="")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
    parser.add_argument("gfa", help="GFA file with graph segment sequences and links")
    parser.add_argument(
        "paths",
        help="Rukki paths file with a GAF-style (>node<node) path in column two",
    )
    parser.add_argument(
        "result",
        nargs="?",
        default="-",
        help="Output FASTA file, or stdout when omitted or '-'",
    )
    parser.add_argument(
        "--trim-flanks-to",
        type=int,
        default=-1,
        help="Trim flanking segments to specified value (disabled if negative)",
    )
    return parser.parse_args(argv)


def main(argv=None):
    from Bio import SeqIO

    args = parse_args(argv)
    if args.trim_flanks_to > 0:
        print(
            "Will trim flanking segments to %dbp" % args.trim_flanks_to,
            file=sys.stderr,
        )

    records = join_paths(args.gfa, args.paths, args.trim_flanks_to)
    if args.result == "-":
        SeqIO.write(records, sys.stdout, "fasta")
    else:
        with open(args.result, "w") as output:
            SeqIO.write(records, output, "fasta")


if __name__ == "__main__":
    main()
