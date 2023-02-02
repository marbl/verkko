#!/usr/bin/env perl

use strict;

my $assembly   = "";
my %contigLength;
my $graph      = "";
my $graphmap   = "";
my $minLength  = 100000;
my %contaminantSeq;

my $threads    = 4;

my $output;

while (scalar(@ARGV) > 0) {
    my $opt = shift @ARGV;

    if    ($opt eq "--assembly") {
        $assembly = shift @ARGV;
    }
    elsif ($opt eq "--graph") {
        $graph = shift @ARGV;
    }
    elsif ($opt eq "--graphmap") {
        $graphmap = shift @ARGV;
    }
    elsif ($opt eq "--minlength") {
        $minLength = int(shift @ARGV);
    }
    elsif ($opt eq "--contaminant") {
        my $prefix   = shift @ARGV;
        my $sequence = shift @ARGV;

        $contaminantSeq{$prefix} = $sequence;
    }
    elsif ($opt eq "--threads") {
        $threads = int(shift @ARGV);
    }
    elsif ($opt eq "--output") {
        $output = shift @ARGV;
    }
    else {
        die "Unknown option '$opt'.\n";
    }
}

if (($assembly eq "") ||
    ($graph    eq "") ||
    (scalar(keys %contaminantSeq) eq 0)) {
    die "usage: $0 ...\n";
}



#  Read sequence lengths into global %contigLength.
#
sub loadSequenceLengths ($) {
    my $a = shift @_;
    my $n = 0;
    my $t = 0;

    print STDERR "Loading sequence lengths from '$a'.\n";

    open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
    while (!eof(F)) {
        my $h = <F>;
        my $s = <F>;

        $h =~ s/^\s+//;   $s =~ s/^\s+//;
        $h =~ s/\s+$//;   $s =~ s/\s+$//;

        if ($h =~ m/^>(\S+)\s*/) {
            $contigLength{$1}  = length($s);
            $n                += 1;
            $t                += length($s);
        }
        else {
            die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
        }
    }
    close(F);

    print STDERR "  Found $n sequences with total length $t bp.\n";
    print STDERR "\n";
}


#  Map in 6-layoutContigs/unitig-popped.layout.scfmap
#  path haplotype1-0000002 mat_from_utig4-117


#  Find disconnected nodes in the graph.  Returns an array of node names.
#
#  cat ../{input.graph} | grep ^L |                 awk '{{ print \$2"\n"\$4 }}' | sort | uniq > ./nodes-with-links
#  cat ../{input.graph} | grep ^S | sed s/LN:i:// | awk '{{ print \$2        }}' | sort | uniq > ./nodes
#
#  diff --side-by-side --suppress-common-lines ./nodes ./nodes-with-links | awk '{{ print \$1 }}' > ./nodes-isolated
#
sub findDisconnected ($$$) {
    my $g = shift @_;
    my $m = shift @_;
    my $minl = shift @_;
    my %gtoc;
    my %nodes;
    my %edges;
    my @discShort;
    my @discLong;

    print STDERR "Loading graph to contig name map from '$m'.\n";

    open(M, "< $m") or die "Failed to open graph-to-contig name map '$m' for reading: $!\n";
    while (<M>) {
        s/^\s+//;
        s/\s+$//;

        if (m/^path\s+(\S+)\s+(\S+_from_|\S+_unused_){0,1}(\S+)$/) {
            my $cn = $1;
            my $gn = $3;
            $gtoc{$gn} = $cn;
        } elsif (m/path/) {
            die "Failed to parse path '$_'\n";
        }
    }
    close(M);

    print STDERR "Finding disconnected nodes in '$g' shorter than $minl.\n";

    open(G, "< $g") or die "Failed to open graph '$g' for reading: $!\n";
    while (<G>) {
        my ($t, $n1, undef, $n2) = split '\t', $_;

        #next   if (!exists($gtoc{$n1}));   #  The graph seems to have stuff not in the output?
        #next   if (!exists($gtoc{$n2}));

        #  Not sure this is correct.  On a rukki run, the warnings below
        #  trigger for some nodes.  I suspect these are the things rukki is
        #  building scaffolds with.

        if    ($t eq "S") {
            my $c1 = $gtoc{$n1};

            #print STDERR "WARNING1: graph name '$n1' not found.\n"   if (!defined($c1));

            $nodes{$c1}++;
        }
        elsif ($t eq "L") {
            my $c1 = $gtoc{$n1};
            my $c2 = $gtoc{$n2};

            #print STDERR "WARNING2: graph name '$n1' not found.\n"   if (!defined($c1));
            #print STDERR "WARNING2: graph name '$n2' not found.\n"   if (!defined($c2));

            $edges{$c1}++;
            $edges{$c2}++;
        }
    }
    close(G);

    foreach my $c (keys %nodes) {
        next   if (! exists($edges{$c}));

        if ($contigLength{$c} <  $minl) {
            push @discShort, $c;
        } else {
            push @discLong, $c;
        }
    }

    printf STDERR " Scanned %7d nodes.\n", scalar(keys %nodes);
    printf STDERR " Scanned %7d edges.\n", scalar(keys %edges);
    printf STDERR " Found %5d short disconnected nodes.\n", scalar(@discShort);
    printf STDERR " Found %5d long  disconnected nodes.\n", scalar(@discLong);
    printf STDERR "\n";

    return @discShort;
}


sub mashMap ($$$$) {
    my $a      = shift @_;
    my $output = shift @_;
    my $cp     = shift @_;
    my $cs     = shift @_;

    my $map;

    if (-e "$output.$cp.mashmap.out") {
        print STDERR "Using pre-computed mashmap results for '$cp' ('$cs').\n";
    }
    else {
        print STDERR "Running mashmap against '$cp' ('$cs').\n";

        $map  = "mashmap";
        $map .= " --ref '$a'";
        $map .= " --query '$cs'";
        $map .= " --perc_identity 95";
        $map .= " --segLength 10000";
        $map .= " --filter_mode none";
        $map .= " --threads $threads";
        $map .= " --output '$output.$cp.mashmap.out'";
        $map .= " > $output.$cp.mashmap.err 2>&1";

        system($map);
    }
}


#
#  Filters mashmap output to remove low-quality and spurious matches.
#  Remembers the best match.
#
#  MashMap output columns:
#     1 query
#     2 len
#     3 bgn
#     4 end
#     5 strand
#     6 ref
#     7 len
#     8 bgn
#     9 end
#    10 ident
#
#  awk '{{if (\$NF > 99 && ((\$9-\$8)/\$7 > 0.50 || (\$9-\$8)/\$7 > 0.25 && \$7 < 50000)) print \$6}}' | sort | uniq > contaminant.list
#
sub filterMash ($) {
    my  $cp = shift @_;

    my $nHits;
    my $nSaveI;
    my $nSave1;
    my $nSave2;

    my  %hits;
    my ($best, $bestI, $bestC) = ("", 0, 0);

    open(F, "< $output.$cp.mashmap.out") or die "Failed to open mashmap output '$output.$cp.mashmap.out' for reading: $!\n";
    while (<F>) {
        my ($q, $qlen, $qbgn, $qend, $s, $ref, $rlen, $rbgn, $rend, $ident) = split '\s+', $_;

        my $rcov = ($rend - $rbgn) / $rlen;

        my $gi = ($ident >= 99.0);
        my $g1 = ($rcov >= 0.50);
        my $g2 = ($rcov >= 0.25) && ($rlen < 50000);

        $nHits++;
        $nSaveI += 1   if ($gi);
        $nSave1 += 1   if ($g1) && (!$gi);
        $nSave2 += 1   if ($g2) && (!$gi);

        if (($gi) && ($g1 || $g2)) {   #  If ident is good and the contig is mostly
            $hits{$ref}++;             #  covered, save the contig as a contaminant.

            if ($ident > $bestI) {
                $best  = $ref;         #  Save this contig if it is the best we've seen.
                $bestI = $ident;
                $bestC = $rcov;
            }
        }
    }
    close(F);

    printf "  Found %5d       hits.\n", $nHits;
    printf "  Found %5d       hits with good identity.\n", $nSaveI;
    printf "  Found %5d       hits with good coverage but bad identity.\n", $nSave1;
    printf "  Found %5d short hits with good coverage but bad identity.\n", $nSave2;
    printf "  Found %5d       contaminant contigs.\n", scalar(keys %hits);
    printf "  Best hit is %6.2f%% identity and %6.2f%% coverage ('$best').\n", $bestI, 100.0 * $bestC;
    printf "\n";

    return($best, keys %hits);
}



#  Extract '@filter' sequences from $a.
#    If $savehits == 1, save the sequences     in '@filter' to the output.
#    If $savehits == 0, save the sequences NOT in '@filter' to the output.
#
sub filterSequences ($$$$$@) {
    my $a           = shift @_;
    my $savehits    = shift @_;
    my $cp          = shift @_;
    my $output      = shift @_;
    my $exemplar    = shift @_;
    my @filter      =       @_;
    my %filter;

    foreach my $f (@filter) {
        $filter{$f}++;
    }

    print STDERR "Filtering '$cp' from '$a'.\n";

    open(O, "> $output.$cp.fasta")          or die "Failed to open '$output.$cp.fasta' for output: $!\n";
    open(E, "> $output.$cp.exemplar.fasta") or die "Failed to open '$output.$cp.exemplar.fasta' for output: $!\n";

    open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
    while (!eof(F)) {
        my $h = <F>;
        my $s = <F>;

        $h =~ s/^\s+//;   $s =~ s/^\s+//;
        $h =~ s/\s+$//;   $s =~ s/\s+$//;

        if ($h =~ m/^>(\S+)\s*/) {
            my $n = $1;
            my $e = exists($filter{$n});

            print O "$h\n$s\n"  if ( $savehits &&  $e);
            print O "$h\n$s\n"  if (!$savehits && !$e);

            print E "$h\n$s\n"  if ($n eq $exemplar);
        }
        else {
            die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
        }
    }
    close(F);

    close(E);
    close(O);

    unlink "$output.$cp.exemplar.fasta" if (!defined($exemplar));

    print STDERR "\n";
}






#
#  Main
#

my %contigLengths = loadSequenceLengths($assembly);
my @disconnected  = findDisconnected($graph, $graphmap, $minLength);

filterSequences($assembly, 1, "disconnected", $output, undef, @disconnected);

my @crud = @disconnected;

foreach my $cp (keys %contaminantSeq) {
    mashMap($assembly, $output, $cp, $contaminantSeq{$cp});

    my ($exemplar, @contaminants) = filterMash($cp);

    filterSequences($assembly, 1, $cp, $output, $exemplar, @contaminants);

    push @crud, @contaminants;
}


filterSequences($assembly, 0, "filtered", $output, undef, @crud);
