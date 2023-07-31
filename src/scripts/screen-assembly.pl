#!/usr/bin/env perl

use strict;

my $assembly   = "";       #  Input file: the contigs after consensus; 7-consensus/unitig-popped.fasta
my %contigLength;          #  Map from contig name to actual length
my $graph      = "";       #  Input file: the graph; 6-rukki/unitig-popped-unitig-normal-connected-tip.noseq.gfa
my $graphmap   = "";       #  Input file: contig layout map; 6-layoutContigs/unitig-popped.layout.scfmap
my %gtoc;                  #  Map from graph name to contig name (from the above file)
my %ctog;                  #  Map from contig name to graph name.
my $hificov    = "";       #  Input file: hifi coverage; 5-untip/unitig-popped-unitig-normal-connected-tip.hifi-coverage.csv
my %hifiCoverage;          #  Map from contig name to coverage (from above file)

my %contaminantSeq;        #  Input contaminant files: map from 'prefix' to 'sequence file'.

my $minLength  = 100000;   #  Parameter: cutoff for being a 'short' contig.
my $threads    = 4;        #  Parameter: number of threads to use when mapping.
my $mashmap    = "mashmap";

my $output;                #  Parameter: output file name prefix

my $fastmode   = 0;


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
    elsif ($opt eq "--hifi-coverage") {
        $hificov = shift @ARGV;
    }
    elsif ($opt eq "--minlength") {
        $minLength = int(shift @ARGV);
    }
    elsif ($opt eq "--contaminant") {
        my $prefix   = $ARGV[0];
        my $sequence = $ARGV[1];

        while (($prefix !~ m/^--/) && (-e $sequence)) {
            $contaminantSeq{$prefix} = $sequence;

            shift @ARGV;
            shift @ARGV;

            $prefix   = $ARGV[0];
            $sequence = $ARGV[1];
        }
    }
    elsif ($opt eq "--threads") {
        $threads = int(shift @ARGV);
    }
    elsif ($opt eq "--output") {
        $output = shift @ARGV;
    }
    elsif ($opt eq "--mashmap") {
        $mashmap = shift @ARGV;
    }
    elsif ($opt eq "--fast") {    #  For debugging the contaminant filtering;
        $fastmode =    1;         #  doesn't write filtered sequence outputs.
    }
    else {
        die "Unknown option '$opt'.\n";
    }
}

if (($assembly eq "") ||
    ($graph    eq "") ||
    ($output   eq "")) {
    print "usage: $0 --assembly X.fasta --graph X.gfa ...\n";
    print "  --assembly X.fasta           sequences to screen for crud\n";
    print "  --graph X.gfa                graph to decide if a short contig is disconnected\n";
    print "  --graphmap X.scfmap          layoutContigs scfmap to rename graph nodes to contig names\n";
    print "  --hifi-coverage X.csv        contig read coverage\n";
    print "  --minlength L                contigs shorter than L bp are deemed 'short' (default: 100000)\n";
    print "  --contaminant N F [N F ...]  label N and contaminant fasta F to screen; example:\n";
    print "                                  --contaminant ebv  ebv.fasta.gz \\\n";
    print "                                                mito mito.fasta   \\\n";
    print "                                                rdna r.fasta.gz\n";
    print "  --threads T                  number of mashmap compute threads to use\n";
    print "  --fast                       run faster; DOES NOT CREATE FASTA OUTPUTS!\n";
    print "  --output X                   write outputs to prefix X\n";
    print "                                  X.disconnected.fasta\n";
    print "                                  X.ebv.exemplar.fasta\n";
    print "                                  X.ebv.fasta\n";
    print "                                  X.ebv.mashmap.err\n";
    print "                                  X.ebv.mashmap.out\n";
    print "                                  X.fasta\n";
    print "                                  X.mito.exemplar.fasta\n";
    print "                                  X.mito.fasta\n";
    print "                                  X.mito.mashmap.err\n";
    print "                                  X.mito.mashmap.out\n";
    print "                                  X.rdna.exemplar.fasta\n";
    print "                                  X.rdna.fasta\n";
    print "                                  X.rdna.mashmap.err\n";
    print "                                  X.rdna.mashmap.out\n";
    exit(1);
}





#  Read sequence lengths into global %contigLength.
#
sub loadSequenceLengths ($) {
    my $a = shift @_;
    my $n = 0;
    my $t = 0;

    print "Loading sequence lengths from '$a'.\n";

    if (-s "$a.lengths") {
        open(F, "< $a.lengths");
        while (<F>) {
            my ($n, $l) = split '\s+', $_;
            $contigLength{$n} = $l;
            $n                += 1;
            $t                += $l;
        }
        close(F);
    }

    else {
        open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
        open(L, "> $a.lengths");

        while (!eof(F)) {
            my $h = <F>;
            my $s = <F>;

            $h =~ s/^\s+//;   $s =~ s/^\s+//;
            $h =~ s/\s+$//;   $s =~ s/\s+$//;

            if ($h =~ m/^>(\S+)\s*/) {
                my $ctg = $1;
                my $len = length($s);

                $contigLength{$ctg}  = $len;
                $n                  +=  1;
                $t                  += $len;

                print L "$ctg\t$len\n"
            }
            else {
                die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
            }
        }
        close(L);
        close(F);
    }

    print "  Found $n sequences with total length $t bp.\n";
    print "\n";
}



#  Read hifi coverage for the contigs.
#
sub loadHifiCoverage ($) {
    my $c = shift @_;

    print "Loading contig coverage from '$c'.\n";

    open(F, "< $c") or die "Failed to open '$c' for reading: $!.\n";
    while (<F>) {
        my ($ident, $cov) = split '\s+', $_;

        if (exists($gtoc{$ident})) {
            #print "$ident -> $gtoc{$ident} -> $cov\n";

            $hifiCoverage{$gtoc{$ident}} = $cov;
        }
        else {
        }
    }
    close(F);

    print "  Found ", scalar(keys %hifiCoverage), " contig coverage values.\n";
    print "\n";
}



#  Load a map from graph-contig-names to assembly-contig-names.
#
sub loadGtoC ($) {
    my $m = shift @_;

    print "Loading graph to contig name map from '$m'.\n";

    open(M, "< $m") or die "Failed to open graph-to-contig name map '$m' for reading: $!\n";
    while (<M>) {
        s/^\s+//;
        s/\s+$//;

        if (m/^path\s+(\S+)\s+(\S+_from_|\S+_unused_){0,1}(\S+)$/) {
            my $cn = $1;
            my $gn = $3;
            $gtoc{$gn} = $cn;
            $ctog{$cn} = $gn;
        } elsif (m/path/) {
            die "Failed to parse path '$_'\n";
        }
    }
    close(M);

    print "  Found ", scalar(keys %gtoc), " contig names.\n";
    print "\n";
}



#  Find disconnected nodes in the graph.  Returns an array of node names.
#
#  cat ../{input.graph} | grep ^L |                 awk '{{ print \$2"\n"\$4 }}' | sort | uniq > ./nodes-with-links
#  cat ../{input.graph} | grep ^S | sed s/LN:i:// | awk '{{ print \$2        }}' | sort | uniq > ./nodes
#
#  diff --side-by-side --suppress-common-lines ./nodes ./nodes-with-links | awk '{{ print \$1 }}' > ./nodes-isolated
#
sub findDisconnected ($$) {
    my $g = shift @_;
    my $minl = shift @_;
    my %nodes;
    my %edges;
    my @discShort;
    my @discLong;

    print "Finding disconnected nodes in '$g' shorter than $minl.\n";

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

            #print "WARNING1: graph name '$n1' not found.\n"   if (!defined($c1));

            $nodes{$c1}++;
        }
        elsif ($t eq "L") {
            my $c1 = $gtoc{$n1};
            my $c2 = $gtoc{$n2};

            #print "WARNING2: graph name '$n1' not found.\n"   if (!defined($c1));
            #print "WARNING2: graph name '$n2' not found.\n"   if (!defined($c2));

            $edges{$c1}++;
            $edges{$c2}++;
        }
    }
    close(G);

    foreach my $c (keys %nodes) {
        next   if (exists($edges{$c}));

        if ($contigLength{$c} <  $minl) {
            push @discShort, $c;
        } else {
            push @discLong, $c;
        }
    }

    printf " Scanned %7d nodes.\n", scalar(keys %nodes);
    printf " Scanned %7d edges.\n", scalar(keys %edges);
    printf " Found %5d short disconnected nodes.\n", scalar(@discShort);
    printf " Found %5d long  disconnected nodes.\n", scalar(@discLong);
    printf "\n";

    return(@discShort);
}



sub mashMap ($$$$) {
    my $a      = shift @_;
    my $output = shift @_;
    my $cp     = shift @_;
    my $cs     = shift @_;

    my $map;

    if (-s "$output.$cp.mashmap.out") {
        print "Using pre-computed mashmap results for '$cp' ('$cs').\n";
    }
    else {
        print "Running mashmap against '$cp' ('$cs').\n";

        $map  = "$mashmap";
        $map .= " --ref '$a'";
        $map .= " --query '$cs'";
        $map .= " --perc_identity 95";
        $map .= " --segLength 5000";
        $map .= " --filter_mode none";
        $map .= " --threads $threads";
        $map .= " --output '$output.$cp.mashmap.out'";
        $map .= " > $output.$cp.mashmap.err 2>&1";

        print "\n";
        print "$map\n";
        print "\n";

        system($map);
    }
}



#
#  Filters mashmap output to remove low-quality and spurious matches.
#  Remembers the best match.
#
#    pick node that contains most reference
#    and has highest depth
#      pick the five/ten highest coverage contigs
#      then pick the highest depth
#
#
#  MashMap output columns:
#     1 contaminant contig
#     2 len
#     3 bgn
#     4 end
#     5 strand
#     6 assembled contig
#     7 len
#     8 bgn
#     9 end
#    10 ident
#
#  awk '{{if (\$NF > 99 && ((\$9-\$8)/\$7 > 0.50 || (\$9-\$8)/\$7 > 0.25 && \$7 < 50000)) print \$6}}' | sort | uniq > contaminant.list
#
sub filterMash ($$) {
    my $cp = shift @_;
    my $g  = shift @_;

    #
    #  Pass 1: find good alignments.
    #
    #    KY962518.1  44838  0     44837 - haplotype1-0000052  80185  4701  50137  99.8464
    #    KY962518.1  44838  10000 44837 - haplotype1-0000052  80185  49651 79582  98.8657
    #
    #  We'll save a hit if it is > 98% identity and covering a significant
    #  chunk of the assembled contig.  With a 10 Kbp minimum match out of
    #  mashmap, this will ignore:
    #    spurious hits in contigs larger than 100 Kbp
    #    full length hits to contigs with more than 10 copies of the contaminant
    #

    my %rawhits;
    my $maxConLen = 0;   #  Max length of a contaminant sequence.

    open(F, "< $output.$cp.mashmap.out") or die "Failed to open mashmap output '$output.$cp.mashmap.out' for reading: $!\n";
    while (<F>) {
        my @s = split '\s+', $_;
        my $l = scalar(@s);
        my ($con, $conlen, $conbgn, $conend, $s, $ctg, $ctglen, $ctgbgn, $ctgend) = @s;
        my $ident = ($l > 10 ? $s[12] : pop(@s));
        $ident =~ s/id:f://g;
        $ident *= 100 if $ident < 1;

        $maxConLen = $conlen   if ($maxConLen < $conlen);

        my $concov = int(0.5 + 100.0 * ($conend - $conbgn) / $conlen);
        my $ctgcov = int(0.5 + 100.0 * ($ctgend - $ctgbgn) / $ctglen);

        my $gi   = ($ident  >= 97.5);   #  Good identity.
        my $g1   = ($ctgcov >= 10);     #  Covers a goodly chunk of the contig.
        my $g2   = ($concov >= 25);     #  Covers a goodly chunk of the contaminant.

        #print "$ctgbgn $ctgend $ctglen $con $conbgn $conend $conlen $ident\n";

        next  if (!$gi);   #  Ignore bad identity.
        next  if (!$g1);   #  Ignore bad coverage.

        #$rawhits{$ctg} = ()   if (!exists($rawhits{$ctg}));

        #print "$ctgbgn $ctgend $ctglen $con $conbgn $conend $conlen $ident - SAVED\n";
        push @{$rawhits{$ctg}}, "$ctgbgn\0$ctgend\0$ctglen\0$con\0$conbgn\0$conend\0$conlen\0$ident";
    }

    #
    #  Pass 2: compute how much of the contig is covered by hits, and if a
    #  large fraction of the contig is covered, flag it as a contaminant.
    #  Then remember the contig with the most read coverage to report as the
    #  exemplar.
    #

    my %hits;
    my $exemplar        = "";
    my $exemplarDepth   = 0;
    my $exemplarBreadth = 0.;

    foreach my $ctg (keys %rawhits) {
        my @rh = sort { $a <=> $b } @{$rawhits{$ctg}};   #  Sort by contig begin position.

        my $covlen = 0;   #  Number of bp we cover in the assembled contig.
        my $covend = 0;   #  Highest position on the contig we've covered.
        my $breadth = 0.;

        foreach my $h (@rh) {
            my ($ctgbgn, $ctgend, $ctglen, $con, $conbgn, $conend, $conlen, $ident) = split '\0', $h;

            if ($covend <= $ctgbgn) {
                $covlen += $ctgend - $ctgbgn;   #  No overlap with an existing hit.
                $covend  = $ctgend;
            }
            elsif ($ctgend <= $covend) {        #  Completely contained in existing hits.
                #covlen += 0;
                #covend  = $covend;
            }
            else {                              #  An extension to what we've covered already.
                $covlen += $ctgend - $covend;
                $covend  = $ctgend;
            }
            $breadth = $covlen / $conlen;
        }
        next   if ($covlen < 0.5 * $contigLength{$ctg});    #  Not covering enough of the contig.

        $hits{$ctg}++;

        #printf "Checking contig %s with dreadth %f and coverage %f and current best is %f and %f\n", $ctg, $breadth,  $hifiCoverage{$ctg}, $exemplarBreadth, $exemplarDepth;
        if (($exemplarBreadth < 0.90 && $breadth > $exemplarBreadth) || ($breadth > 0.90 && $hifiCoverage{$ctg} > $exemplarDepth) || ($breadth == $exemplarBreadth &&  $hifiCoverage{$ctg} > $exemplarDepth)) {
            #  lower breadth of coverage of the reference than the existing exemplar
            #  or covers the full reference and has more coverage
            #  or same breadth but lower depth

            $exemplar        = $ctg;
            $exemplarDepth   = $hifiCoverage{$ctg};
            $exemplarBreadth = $breadth;
        }
    }

    sub sumLen (@) {
        my $sum = 0;
        foreach (@_) {
            $sum += $contigLength{$_};
        }
        return $sum;
    }

    #printf "  Found %5d       hits.\n", $nHits;
    #printf "  Found %5d       hits with good identity.\n", $nSaveI;
    #printf "  Found %5d       hits with good coverage but bad identity.\n", $nSave1;
    #printf "  Found %5d short hits with good coverage but bad identity.\n", $nSave2;
    printf "  Found %5d contaminant contigs of total length %d.\n", scalar(keys %hits), sumLen(keys %hits);


    #  Load graph edges, ignoring orientation.
    #  See comments in findDisconnected() about missing contig name warnings.

    my %edges;
    my %lost;

    open(G, "< $g") or die "Failed to open graph '$g' for reading: $!\n";
    while (<G>) {
        my ($t, $n1, undef, $n2) = split '\t', $_;

        if ($t eq "L") {
            my $c1 = $gtoc{$n1};
            my $c2 = $gtoc{$n2};

            if (defined($c1) && defined($c2)) {
                push @{$edges{$c1}}, $c2;
                push @{$edges{$c2}}, $c1;
            } else {
                $lost{$n1} = 1   if (!defined($c1));
                $lost{$n2} = 1   if (!defined($c2));
            }
        }
    }
    close(G);

    #  get_layout_from_mbg.py discards graph contigs from consensus if they
    #  have no reads assigned.  This will warn about those.
    #
    #foreach my $c (sort keys %lost) {
    #    print "No contig found for '$c'.\n";
    #}

    #  Reset %hits to make it represent the number of hops from an original contig hit.

    foreach my $ctg (keys %hits) {
        $hits{$ctg} = 0;
    }

    #  Iterate over the hits, adding anything adjacent (up to four hops from
    #  an origianlly identified contaminant contig) if it is of comparable
    #  size.

    for (my $d=1; $d <= 4; $d++) {
        my %extended;

        print "    extension pass $d:\n";

        foreach my $ctg (keys %hits) {
            foreach my $adj (@{$edges{$ctg}}) {
                next   if (exists($hits{$adj}));                    #  Skip if we're already in the set.
                next   if ($contigLength{$adj} > 4 * $maxConLen);   #  Skip if the adj contig is much bigger.
                next   if (exists($extended{$adj}));                #  Skip if we've already got it.

                printf "      %20s / %-20s --> %20s / %s\n", $ctg, $ctog{$ctg}, $adj, $ctog{$adj};

                $extended{$adj} = $hits{$ctg} + 1;
            }
        }

        foreach my $ctg (keys %extended) {
            if (!exists($hits{$ctg})) {
                $hits{$ctg} = $extended{$ctg};
            }
        }
    }

    printf "  Found %5d contaminant contigs of total length %d (after extending).\n", scalar(keys %hits), sumLen(keys %hits);
    printf "  Best hit has %6.2fx read coverage. ('%s').\n", $exemplarDepth, $exemplar;
    printf "\n";

    return($exemplar, keys %hits);
}




#  Extract '@filter' sequences from $a.
#
#   - Sequences in @filter are copied to the output file.
#   - If $exemplar is defined, that specific sequence is saved, too.
#
sub filterSequences ($$$$@) {
    my $a           = shift @_;
    my $output      = shift @_;
    my $cp          = shift @_;
    my $exemplar    = shift @_;
    my @filter      =       @_;
    my %filter;

    #  Make lookups of sequences we want to filter easier.
    open(O, "> $output.ids") or die "Failed to open '$output.ids' for output: $!\n";
    foreach my $f (@filter) {
        print O "$f\t$ctog{$f}\n";
        $filter{$f}++;
    }
    close(O);

    return   if ($fastmode);

    #  Open some output files.
    if (defined($exemplar)) {
        print "Filtering '$cp' from '$a'.\n";
        open(O, "> $output.fasta")          or die "Failed to open '$output.fasta' for output: $!\n";
        open(E, "> $output.exemplar.fasta") or die "Failed to open '$output.exemplar.fasta' for output: $!\n";
    }
    else {
        print "Extracting sequences from '$a'.\n";
        open(O, "> $output.fasta")          or die "Failed to open '$output.fasta' for output: $!\n";
    }

    #  Scan the input, copying desired sequences to the output.
    open(F, "< $a") or die "Failed to open '$a' for reading: $!\n";
    while (!eof(F)) {
        my $h = <F>;
        my $s = <F>;

        $h =~ s/^\s+//;   $s =~ s/^\s+//;
        $h =~ s/\s+$//;   $s =~ s/\s+$//;

        if ($h =~ m/^>(\S+)\s*/) {
            my $n = $1;

            print O "$h\n$s\n"  if (exists($filter{$n}));
            print E "$h\n$s\n"  if ($n eq $exemplar) && (defined($exemplar));
        }
        else {
            die "Failed to find ident line in input line '", substr($h, 0, 40), "'.\n";
        }
    }
    close(F);

    close(E)   if (defined($exemplar));
    close(O);

    unlink "$output.exemplar.fasta" if (!defined($exemplar));

    print "\n";
}






#
#  Main
#

loadSequenceLengths($assembly);
loadGtoC($graphmap);

my @disconnected  = findDisconnected($graph, $minLength);
filterSequences($assembly, "$output.disconnected", "disconnected", undef, @disconnected);

my @crud = @disconnected;

loadHifiCoverage($hificov);

foreach my $cp (sort keys %contaminantSeq) {
    mashMap($assembly, $output, $cp, $contaminantSeq{$cp});

    my ($exemplar, @contaminants) = filterMash($cp, $graph);

    filterSequences($assembly, "$output.$cp", $cp, $exemplar, @contaminants);
    push @crud, @contaminants;
}

my %gold = %contigLength;
my @gold;

foreach my $k (@crud) {         #  Remove the crud from the gold.
    delete $gold{$k};
}
foreach my $k (keys %gold) {    #  Copy gold to a list.
    push @gold, $k;
}

filterSequences($assembly, $output, undef, undef, @gold);

exit(0);
